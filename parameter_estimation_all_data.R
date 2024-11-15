# 必要なライブラリの読み込み
library(deSolve)
library(ggplot2)

# データの読み込み


data = read.csv("Japanese_mpox_data_week_trimmed_SM.csv")  # データをdata.csvから読み込む

# SEIRモデルの定義
N = 125050000
SEIR_model = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = - S * (beta*(I/N))
    dE = beta*(I/N) * S - E * (v + alpha)
    dI = v * E - I * (delta + rho)
    dR = rho * I + alpha * E
    return(list(c(dS, dE, dI, dR)))
  })
}

# 初期パラメータと初期状態の設定
# initial_parameters = c(mu=3.294617e-01, delta=0, rho=6.650022e-01, v=1, beta=0.720308793132039, alpha=0)
initial_parameters = c(mu=0.001, delta=0.821436749390314, rho=0.698359554114671, v=0.776335339245849, beta=0.785748427237125, alpha=0.197112080032751)
# initial_parameters = c(mu=0, delta=0.1, rho=1, v=0.616923966929326, beta=0.604797391016171, alpha=0.2)
initial_state = c(S=54709945, E=50, I=5, R=0)

# 最小二乗法を使用してパラメータ推定を行う関数
fit_function = function(parameters) {
  # モデルの実行
  mod = ode(y=initial_state, times=data$time_month, func=SEIR_model, parms=parameters)

  # モデルの出力と観測データとの誤差を計算
  model_output = as.data.frame(mod)
  error = sum((model_output[, "I"] - data$cases)^2)  # I（感染者数）の誤差を最小化

  return(error)
}

# 最小二乗法によるパラメータ推定
fit = optim(par = initial_parameters, fn = fit_function, method = "L-BFGS-B", lower = c(Lambda=0, delta=0, rho=0.8, v=0, beta=0, alpha=0), upper = c(Lambda=2, delta=0.1, rho=1, v=1, beta=1, alpha=0.2))
# fit = optim(par = initial_parameters, fn = fit_function, method = "L-BFGS-B")

# 推定されたパラメータを表示
estimated_parameters = fit$par
print(estimated_parameters)
print("推定されたパラメータ:")
print(paste0("delta=", estimated_parameters["delta"], ", rho=", estimated_parameters["rho"], ", nu=", estimated_parameters["v"], ", beta=", estimated_parameters["beta"], ", alpha=", estimated_parameters["alpha"], ", lambda=", estimated_parameters["lambda"]))
