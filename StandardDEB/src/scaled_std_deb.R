scaled_std_deb <- function(t, states, params){
  
  L_m = params["L_m"]
  E_G = params["E_G"]
  v = params["v"]
  kap = params["kap"]
  k_J = params["k_J"]
  E_Hp = params["E_Hp"]
  p_Am = params["p_Am"]
  f = params["f"]
  
  E_m = p_Am/ v
  g = E_G/ kap/ E_m;
  k_M = v / (g * L_m)
  l_T = 0
  k = k_J / k_M
  V_m =(v/(k_M*g))^3
  
  u_Hp = E_Hp/(g*E_m*V_m)
  
  
  #derivatives  
  de = k_M * g * (f - states[1]) / states[1]
  dl = k_M / 3 * (states[1] - states[2] - l_T) / (1 + states[1] / g)
  duH = k_M * (states[3] < u_Hp) * ((1 - kap) * states[1] * states[2]^2 * (g + states[2] + l_T) / (g + states[1]) - k * states[3])
  duR = k_M * (states[3] > u_Hp) * ((1 - kap) * states[1] * states[2]^2 * (g + states[2] + l_T) / (g + states[1]) - k * u_Hp)
  
  #observable quantities
  
  
  
  
  
  
  #return derivatives
  list(c(de, dl, duH, duR))
}