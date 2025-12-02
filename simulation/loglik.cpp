#include <TMB.hpp>
#include <vector>
#include <cmath>   // for fabs
#include <cppad/cppad.hpp>

// helper to clamp Type
template<class Type>
Type clamp(Type x, Type lo, Type hi){
  if(x < lo) return lo;
  if(x > hi) return hi;
  return x;
}

template<class Type>
bool eq_tol(Type x, Type y, Type tol = 1e-8){
  return CppAD::abs(x - y) < tol;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(t);     // observation times
  DATA_VECTOR(M);     // observed M
  DATA_VECTOR(Y);     // observed Y
  int n = t.size();
  
  
  PARAMETER(p1_raw); PARAMETER(p2_raw); PARAMETER(p4_raw); PARAMETER(r_raw);
  
  PARAMETER(logit_c1); PARAMETER(logit_c2); PARAMETER(logit_c4);
  PARAMETER(logit_m1); PARAMETER(logit_m2); PARAMETER(logit_m4);
  
  Type p1 = invlogit(p1_raw);
  Type p2 = invlogit(p2_raw);
  Type p4 = invlogit(p4_raw);
  Type r  = invlogit(r_raw);
  
  Type Uc = 100; 
  Type Um = 100;
  
  Type c1 = Uc * invlogit(logit_c1);
  Type c2 = Uc * invlogit(logit_c2);
  Type c4 = Uc * invlogit(logit_c4);
  
  Type m1 = Um * invlogit(logit_m1);
  Type m2 = Um * invlogit(logit_m2);
  Type m4 = Um * invlogit(logit_m4);
  

  std::vector<int> hidden_states(1, 0);  // initial hidden state
  vector<Type> alpha_prev(1); alpha_prev(0) = 1.0;
  Type nll = 0.0;
  
  for(int step = 1; step < n; step++){
    
    Type deltaT = t(step) - t(step-1);
    Type deltaY = Y(step) - Y(step-1);
    Type deltaM = M(step) - M(step-1);
    
    // event probabilities, clamped
    Type p1_t = clamp(p1 / (1.0 + c1*pow(t(step) - m1, 2)), Type(1e-12), Type(1.0-1e-12));
    Type p2_t = clamp(p2 / (1.0 + c2*pow(t(step) - m2, 2)), Type(1e-12), Type(1.0-1e-12));
    Type p4_t = clamp(p4 / (1.0 + c4*pow(t(step) - m4, 2)), Type(1e-12), Type(1.0-1e-12));
    Type p3_t = clamp(1.0 - p1_t - p2_t - p4_t, Type(1e-12), Type(1.0-1e-12));
    
    // rates
    int K = hidden_states.size();
    vector<Type> rates(K), exp_term(K);
    for(int i=0; i<K; i++){
      Type rate = r * (M(step-1) - Type(hidden_states[i]));
      if(!CppAD::isfinite(rate) || rate < 1e-8) rate = Type(1e-8);
      rates(i) = rate;
      exp_term(i) = dexp(deltaT, rate, false);
    }
    
    vector<Type> alpha_new;
    
    // Event 2: deltaY==1, deltaM==0
    if(eq_tol(deltaY, Type(1.0)) && eq_tol(deltaM, Type(0.0))){
      alpha_new = alpha_prev.array() * exp_term.array() * p2_t;
    }
    // Event 3: deltaY==2, deltaM==-1
    else if(eq_tol(deltaY, Type(2.0)) && eq_tol(deltaM, Type(-1.0))){
      alpha_new = alpha_prev.array() * exp_term.array() * p3_t;
      
      // prune illegal hidden states
      std::vector<int> valid_idx;
      for(int i=0;i<K;i++){
        if(Type(hidden_states[i]) <= M(step)+1e-8) valid_idx.push_back(i);
      }
      if(valid_idx.size()==0) return Type(1e12);
      
      vector<Type> tmp(valid_idx.size());
      std::vector<int> new_states(valid_idx.size());
      for(size_t i=0;i<valid_idx.size();i++){
        tmp(i) = alpha_new(valid_idx[i]);
        new_states[i] = hidden_states[valid_idx[i]];
      }
      alpha_new = tmp;
      hidden_states = new_states;
    }
    // Event 1 or 4: deltaY==0, deltaM==1
    else if(eq_tol(deltaY, Type(0.0)) && eq_tol(deltaM, Type(1.0))){
      vector<Type> diag_part = alpha_prev.array() * exp_term.array() * p1_t;
      vector<Type> shift_part = alpha_prev.array() * exp_term.array() * p4_t;
      
      int L = diag_part.size();
      if(L==0) return Type(1e12);
      
      alpha_new = vector<Type>(L+1);
      for(int i=0;i<L;i++) alpha_new(i) = diag_part(i);
      for(int i=0;i<L;i++) alpha_new(i+1) += shift_part(i);
      
      // add new hidden state
      std::vector<int> new_states(L+1);
      for(int i=0;i<L;i++) new_states[i] = hidden_states[i];
      new_states[L] = new_states[L-1]+1;
      hidden_states = new_states;
      
      // prune invalid states
      std::vector<int> valid_idx;
      for(size_t i=0;i<hidden_states.size();i++){
        if(Type(hidden_states[i]) <= M(step)+1e-8) valid_idx.push_back(i);
      }
      if(valid_idx.size()==0) return Type(1e12);
      
      vector<Type> tmp(valid_idx.size());
      std::vector<int> states2(valid_idx.size());
      for(size_t i=0;i<valid_idx.size();i++){
        tmp(i) = alpha_new(valid_idx[i]);
        states2[i] = hidden_states[valid_idx[i]];
      }
      alpha_new = tmp;
      hidden_states = states2;
    }
    else return Type(1e12);
    
    Type s = alpha_new.sum();
    if(!CppAD::isfinite(s) || s <= 0) return Type(1e12);
    
    nll -= log(s);
    alpha_prev = alpha_new / s;
  }
  
  return nll;
}
