library(ggplot2)
library(gridExtra)
library(cowplot)
library(deSolve)
library(reshape2)
library(FME)
#delta : death rate of Mpox
#rho : Recovery rate
#nu : Symptom rate
#beta : Infectious rate
#alpha : Through(Skip) Symptom rate
#gamma : Vaccination rate

N = 125050000
SEIR_model = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dS = -S*(beta*(I/N)+mu + gamma_1)
        dE = beta*(I/N)*S - E*(mu+nu+alpha + gamma_2)
        dI = nu*E - I*(mu+delta+rho)
        dR = rho*I + alpha*E - R*(mu)
        dV = gamma_1*S + gamma_2*E
        list(c(dS, dE, dI, dR, dV))
    })
}
SEIR_model_normal = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dS = -S*(beta*(I/N)+mu)
        dE = beta*(I/N)*S - E*(mu+nu+alpha)
        dI = nu*E - I*(mu+delta+rho)
        dR = rho*I + alpha*E - R*(mu)
        list(c(dS, dE, dI, dR))
    })
}

initial_0 = c(S=125050000-50-5-70340000, E=50, I=5, R=0) # E=50の時とE=0の時をここで設定してください。
# 以下のパラメータのいずれかをEの値に合わせてコメント解除して実行してください。
# parameters_0 = c(mu=0, delta=0.1, rho=1, nu=0.616924082661425, beta=0.604797385420912, alpha=0.2) # E>Iの決定版(E=50)
# parameters_0 = c(mu=0, delta=0, rho=0.8, nu=1, beta=1, alpha=0) # E<Iの決定版(E=0)
times_0 = seq(0.25, 4.75, 0.01)
out_0 = ode(y=initial_0, times=times_0, func=SEIR_model_normal, parms=parameters_0)
out_0 = data.frame(out_0)
out_0 = out_0[1:2500,]
print(R0_calc_normal(parameters_0["mu"], parameters_0["delta"], parameters_0["rho"], parameters_0["nu"], parameters_0["beta"], parameters_0["alpha"]))

initial_0.5 = c(S=125050000-50-5-70340000, E=0, I=5, R=0, V=70340000) # E=50の時とE=0の時をここで設定してください。
# 以下のパラメータのいずれかをEの値に合わせてコメント解除して実行してください。
# parameters_0.5 = c(mu=0, delta=0.1, rho=1, nu=0.616924082661425, beta=0.604797385420912, alpha=0.2, gamma_1=0.7, gamma_2=0.1) # E>Iの決定版
# parameters_0.5 = c(mu=0, delta=0, rho=0.8, nu=1, beta=1, alpha=0, gamma_1=0.7, gamma_2=0.1) # E<Iの決定版
times_0.5 = seq(0.25, 4.75, 0.01)
out_0.5 = ode(y=initial_0.5, times=times_0.5, func=SEIR_model, parms=parameters_0.5)
out_0.5 = data.frame(out_0.5)
out_0.5 = out_0.5[1:2500,]
print(R0_calc_vaccination(parameters_0.5["mu"],parameters_0.5["delta"], parameters_0.5["rho"], parameters_0.5["nu"], parameters_0.5["alpha"], parameters_0.5["lambda"], parameters_0.5["gamma"]))

initial_0.7 = c(S=125050000-50-5-70340000, E=0, I=5, R=0, V=70340000) # E=50の時とE=0の時をここで設定してください。
# 以下のパラメータのいずれかをEの値に合わせてコメント解除して実行してください。
# parameters_0.7 = c(mu=0, delta=0.1, rho=1, nu=0.616924082661425, beta=0.604797385420912, alpha=0.2, gamma_1=0.1, gamma_2=0.7) # E>Iの決定版
# parameters_0.7 = c(mu=0, delta=0, rho=0.8, nu=1, beta=1, alpha=0, gamma_1=0.1, gamma_2=0.7) # E<Iの決定版
times_0.7 = seq(0.25, 4.75, 0.01)
out_0.7 = ode(y=initial_0.7, times=times_0.7, func=SEIR_model, parms=parameters_0.7)
out_0.7 = data.frame(out_0.7)
out_0.7 = out_0.7[1:2500,]
print(R0_calc_vaccination(parameters_0.7["mu"],parameters_0.7["delta"], parameters_0.7["rho"], parameters_0.7["nu"], parameters_0.7["alpha"], parameters_0.7["lambda"], parameters_0.7["gamma"]))

df_japan = read.csv("Japanese_mpox_data_week_trimmed_SM.csv")

df_normal = read.csv("Fitting_predicted_data_as_week.csv")

# Iのグラフ表示
ggplot(NULL) +
  geom_bar(data=df_japan, aes(x=df_japan$time_month, y=df_japan$cases), stat="identity", fill="#2297E6", alpha=0.5) +

  geom_line(data=out_0, aes(out_0$time, out_0$I, colour="Normal"), size=1.5) +
  geom_line(data=out_0.5, aes(out_0.5$time, out_0.5$I, colour="10%"), size=1.5) +
  geom_line(data=out_0.7, aes(out_0.7$time, out_0.7$I, colour="70%"), size=1.5) +
  scale_color_manual(name = "Percentage of\n post-exposure\nprophylaxis", values = c("Normal"="#DF536B", "10%"="#FF9900", "70%"="#9A0079")) +
  theme(legend.position="top") +
  theme_classic(base_size=16) +
  ylab("Infected Population") +
  xlab("time(months)")

# Eのグラフ表示
# ggplot(NULL) +
#   # geom_bar(data=df_japan, aes(x=df_japan$time_month, y=df_japan$cases), stat="identity", fill="#2297E6", alpha=0.5) +
# 
#   geom_line(data=out_0, aes(out_0$time, out_0$S, colour="Normal"), size=1.5) +
#   geom_line(data=out_0.5, aes(out_0.5$time, out_0.5$S, colour="10%"), size=1.5) +
#   geom_line(data=out_0.7, aes(out_0.7$time, out_0.7$S, colour="70%"), size=1.5) +
#   scale_color_manual(name = "Percentage of\n post-exposure\nprophylaxis", values = c("Normal"="#DF536B", "10%"="#FF9900", "70%"="#9A0079")) +
#   theme(legend.position="top") +
#   theme_classic(base_size=16) +
#   ylab("Susceptible Population") +
#   xlab("time(months)")

# Iのグラフ表示
# ggplot(NULL) +
#   # geom_bar(data=df_japan, aes(x=df_japan$time_month, y=df_japan$cases), stat="identity", fill="#2297E6", alpha=0.5) +
# 
#   geom_line(data=out_0, aes(out_0$time, out_0$E, colour="Normal"), size=1.5) +
#   geom_line(data=out_0.5, aes(out_0.5$time, out_0.5$E, colour="10%"), size=1.5) +
#   geom_line(data=out_0.7, aes(out_0.7$time, out_0.7$E, colour="70%"), size=1.5) +
#   scale_color_manual(name = "Percentage of\n post-exposure\nprophylaxis", values = c("Normal"="#DF536B", "10%"="#FF9900", "70%"="#9A0079")) +
#   theme(legend.position="top") +
#   theme_classic(base_size=16) +
#   ylab("Exposure Population") +
#   xlab("time(months)")

# Rのグラフ表示
# ggplot(NULL) +
#   # geom_bar(data=df_japan, aes(x=df_japan$time_month, y=df_japan$cases), stat="identity", fill="#2297E6", alpha=0.5) +
# 
#   geom_line(data=out_0, aes(out_0$time, out_0$R, colour="Normal"), size=1.5) +
#   geom_line(data=out_0.5, aes(out_0.5$time, out_0.5$R, colour="10%"), size=1.5) +
#   geom_line(data=out_0.7, aes(out_0.7$time, out_0.7$R, colour="70%"), size=1.5) +
#   scale_color_manual(name = "Percentage of\n post-exposure\nprophylaxis", values = c("Normal"="#DF536B", "10%"="#FF9900", "70%"="#9A0079")) +
#   theme(legend.position="top") +
#   theme_classic(base_size=16) +
#   ylab("Recovered Population") +
#   xlab("time(months)")


