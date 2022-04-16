# install required package names
pkgs = c("dplyr", "shape")
pkgs_uninstalled = pkgs[!( pkgs %in% installed.packages()[,"Package"] )]
if (length(pkgs_uninstalled)) {
  install.packages(pkgs_uninstalled)
}

# include packages
lapply(pkgs, require, character.only = TRUE)


# 最大人數
Wmax = 6
# b的比例
ratio = 0.1
# U型生產線之速度
VC = 1
# 移動速度
VW = 4

#大4-17(17)
relation = matrix(c(1,2,2,3,3,4,4,5,5,6,6,7,6,9,7,8,8,10,9,10,10,11,10,12,11,13,12,13,13,14,13,15,14,16,15,17,16,17,17,18,18,19,19,20,20,21,21,22,22,23,22,24,22,25,23,29,24,28,25,26,25,27,28,29,26,29,27,29), ncol = 2, byrow = T)
T = c(23,21,26,32,13,24,16,14,22,5,16,35,55,5,11,5,39,8,22,34,31,64,61,24,20,15,9,48,31)
      
#大4-18(18)
#relation = matrix(c(1,2,2,3,3,4,4,6,5,6,6,7,7,8,8,9,9,10,10,12,11,12,12,13,13,15,14,15,15,16,15,17,16,18,17,19,18,19,19,20,19,23,20,21,20,22,21,25,22,25,23,24,24,25,25,26,26,27,27,28,27,29,27,30,28,34,29,33,30,31,30,32,31,34,32,34,33,34), ncol = 2, byrow = T)
#T = c(23,21,26,23,11,32,13,24,16,14,22,5,16,35,55,5,11,5,39,12,11,15,8,22,34,31,64,61,24,20,15,9,48,31)
            
#大4-19(19)
#relation = matrix(c(1,2,2,3,3,6,4,6,5,6,6,7,7,8,8,9,8,10,9,11,10,11,11,12,11,13,12,14,13,14,14,15,14,16,15,17,15,18,15,19,16,17,16,18,16,22,17,22,18,22,19,20,20,21,21,22,22,23,22,24,22,26,22,27,23,30,24,25,25,30,26,29,27,28,28,30,29,30,30,31,31,32,32,33,33,34,33,35,33,36,34,40,35,39,36,37,36,38,37,40,38,40,39,40), ncol = 2, byrow = T)
#T = c(23,21,26,45,20,32,31,22,22,61,5,16,35,55,5,11,40,42,55,22,5,39,12,14,20,8,25,36,22,34,82,31,64,61,24,20,15,9,48,31)

#處理先行關係矩陣
PR = matrix(0, nrow = max(relation), ncol = max(relation))
PR[relation] = 1

#計算b
b = sum(T)*ratio

#計算a
a = 0.5*(sum(T)/VC - b)

#控制小數點呈現位數
specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))

#函數：產生初始解
solGen = function(Wmax, PR, T, a, b){
  
  J = length(T)
  j = 1
  X = {}
  T_rank = {}
  exclude = {}
  PR_modi = PR
  while (j < J){
    
    C = which(colSums(PR_modi)== 0)
    #print(C)
    if (j>1){
      C = C[!(C %in% exclude)]
      #print(C)
    }
    if (length(C)==1){
      i = C
    }else{
      i = sample(C,1)
    }
    #print(i)
    
    exclude = c(exclude, i)
    
    X = c(X, i)
    #print(X)
    T_rank = c(T_rank,T[i])
    PR_modi = PR[-X,]
    j = j + 1
  }
  
  i = (1:J)[!((1:J) %in% exclude)]
  X = c(X, i)
  T_rank = c(T_rank,T[i])

  # location -> operator
  Y = rep(0, J)
  j_overbar = X[1]
  j_underbar = X[J]
  D = c(j_overbar,j_underbar)
  k = 1
  k1 = 1
  k2 = J

  while (k <= Wmax){
    totaltime = 0
    while (totaltime < sum(T)/Wmax){
      if (j_overbar != j_underbar){
        i = sample(D,1)
        #print(i)
        Y[i] = k
        #print(Y)
        totaltime = totaltime + T[i]
        #print(totaltime)
        if (i == j_overbar){
          k1 = k1 + 1
          j_overbar = X[k1]
        }else{
          k2 = k2 - 1
          j_underbar = X[k2]
        }
        D = c(j_overbar,j_underbar)
      }else{
        Y[j_overbar] = Wmax
        totaltime = sum(T)/Wmax # 強迫跳出
      }
    }
    k = k + 1
  }
  
  output = {}
  temp = rbind(X, 1:J, T_rank)
  
  for (i in 1:J){
    output = cbind(output, temp[,which(X==i)])
  }
  output = rbind(output, Y)
  
  row.names(output)=c("index", "X", "T", "Y")
  #return(list(output=output, T_rank=T_rank))
  

  # 計算U型的上水平線
  Tx = 0
  i = 1
  col_coordinate = c(0,b)
  while (Tx < a){
    Tx = Tx + output['T',][which(output['X',]==i)]
    i = i + 1
    col_coordinate = rbind(col_coordinate, c(Tx, b))
  } 
  # 將X軸超過者扣掉，再分配給Y軸
  Tx_over = col_coordinate[nrow(col_coordinate),1] - a
  col_coordinate[nrow(col_coordinate),] = col_coordinate[nrow(col_coordinate),] - Tx_over
  
  # 計算U型的右垂直線 
  Ty = col_coordinate[nrow(col_coordinate), 2]
  while (Ty > 0){
    Ty = Ty - output['T',][which(output['X',]==i)]
    i = i + 1
    col_coordinate = rbind(col_coordinate, c(a, Ty))
  } 
  
  Ty_over = -col_coordinate[nrow(col_coordinate),2]
  col_coordinate[nrow(col_coordinate),] = c(col_coordinate[nrow(col_coordinate),1] - Ty_over, col_coordinate[nrow(col_coordinate),2] + Ty_over)
  
  Tx = col_coordinate[nrow(col_coordinate), 1]
  while (Tx > 0 & i <= J){
    Tx = Tx - output['T',][which(output['X',]==i)]
    i = i + 1
    col_coordinate = rbind(col_coordinate, c(Tx, 0))
  } 
  
  output = {}
  temp = rbind(X, 1:J, T_rank, t(col_coordinate)[,1:(ncol(t(col_coordinate))-1)], t(col_coordinate)[,2:ncol(t(col_coordinate))])
  
  for (i in 1:J){
    output = cbind(output, temp[,which(X==i)])
  }
  output = rbind(output, Y)
  
  row.names(output)=c("index", "X", "T", "Startx", "Starty", "Endx", "Endy", "Y")
  output
  
  return(output)
}

#函數：計算目標函數
evaluate = function(output, Wmax, a, b, VW, T){

  # 針對各個工作站計算工作時間與移動時間，並計算該布置的周期時間
  CT = {}
  for (i in 1:Wmax){
    
    output_temp = as.matrix(output[, which(output["Y",]==i)])
    output_order = order(output_temp["X",])
    col_start_end = {}
    
    for (j in 1:ncol(output_temp)){
      start = unname(output_temp[c("Startx","Starty"),output_order[j]])
      end = unname(output_temp[c("Endx","Endy"),output_order[j]])
      # 第一個位置，col_start_end第三欄1代表與上一筆有重複、反之則為0
      if (j == 1){
        col_start_end = unname(rbind(col_start_end, c(start,0)))
        col_start_end = unname(rbind(col_start_end, c(end,1)))
        # 第二個位置之後、且最後一個座標與新的座標重複時
      }else if (j > 1 & sum(col_start_end[nrow(col_start_end),1:2]==start)==2){
        col_start_end = unname(rbind(col_start_end, c(end,1)))
        # 第二個位置之後、且最後一個座標與新的座標無重複時  
      }else if (j > 1 & sum(col_start_end[nrow(col_start_end),1:2]==start)!=2){
        col_start_end = unname(rbind(col_start_end, c(start,0)))
        col_start_end = unname(rbind(col_start_end, c(end,1)))
      }
    }
    # 補充U型生產線的轉折點(a,0)或(a,b)
    col_start_end_mend = {}
    for (k1 in 2:nrow(col_start_end)){
      if (col_start_end[k1-1,1] == a & col_start_end[k1,2] == 0 & col_start_end[k1,3] == 1){
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1-1,])
        col_start_end_mend = rbind(col_start_end_mend, c(a,0,1))
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1,])
      } else if (col_start_end[k1-1,2] == b & col_start_end[k1,1] == a & col_start_end[k1,3] == 1){
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1-1,])
        col_start_end_mend = rbind(col_start_end_mend, c(a,b,1))
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1,])
        # --o--|
        #      |
        # ---o-|        
      } else if (col_start_end[k1-1,1] != a & col_start_end[k1,1] != a & col_start_end[k1-1,2] == b & col_start_end[k1,2] == 0 & col_start_end[k1,3] == 1){
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1-1,])
        col_start_end_mend = rbind(col_start_end_mend, c(a,b,1))
        col_start_end_mend = rbind(col_start_end_mend, c(a,0,1))
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1,])        
      } else {
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1-1,])
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k1,])
      }
    }
    
    #移除重複
    col_start_end_mend = unique(col_start_end_mend)
    #補原點
    col_start_end_mend = rbind(col_start_end_mend, c(unname(output_temp[c("Startx","Starty"),output_order[1]]),0))
    #擴增欄位
    col_start_end_mend_aug = cbind(col_start_end_mend,0,0,0)
    for (k2 in 2:nrow(col_start_end_mend_aug)){
      #若不在U型路線上，計算步行距離/VW=步行時間,存於第五欄
      if (col_start_end_mend_aug[k2,3]==0){
        col_start_end_mend_aug[k2,5]=sqrt(sum((col_start_end_mend_aug[k2-1,1:2]-col_start_end_mend_aug[k2,1:2])^2))/VW
      #若在U型路線上，計算工作距離/VC=工作時間,存於第四欄
      }else{
        col_start_end_mend_aug[k2,4]=sqrt(sum((col_start_end_mend_aug[k2-1,1:2]-col_start_end_mend_aug[k2,1:2])^2))/VC
      }
      #兩者加總存於第五欄
      col_start_end_mend_aug[k2,6]=col_start_end_mend_aug[k2,4]+col_start_end_mend_aug[k2,5]
    }
    #所有時間加總
    CT = c(CT, sum(col_start_end_mend_aug[,6]))
  }
  #1:Wmax中CT最大的直當成績效，值越低越好
  performance = max(CT)
  bottleneck = print(which(CT == max(CT)))
  #print(output)
  #print(diff(output['Y',], lag=1))
  print(performance)
  
  
  return(list(performance=performance, CT = CT, bottleneck=bottleneck))
}






#函數：工作流程的繪圖
plotU = function(output, a, b){
  
  if (mean(output['T',]) <= 1){
    tol = 0.03
  }else if (mean(output['T',]) <= 5 & mean(output['T',]) > 1){
    tol = 0.35
  }else if (mean(output['T',]) <= 10 & mean(output['T',]) > 5){
    tol = 0.65
  }else if (mean(output['T',]) <= 15 & mean(output['T',]) > 10){
    tol = 1
  }else{
    tol = 2
  }

  plot(0, 0, col = "white", xlab = "", ylab = "", xlim = c(0,a+tol), ylim=c(0,b+tol), asp = 1)
  for (ind in 1:ncol(output)){
    i = which(output['X',]==ind)
    slope = (output['Starty',i]-output['Endy',i])/(output['Startx',i]-output['Endx',i])
    if (slope==0 | is.infinite(as.numeric(slope))){
      #
      # 兩點斜率為0或無窮大時
      # -o-o-|   -----|
      #      |        o
      # -----|   -----o
      #
      segments(output['Startx',i], output['Starty',i], output['Endx',i], output['Endy',i], col=output['Y',i], lwd = 10)
      text(x=(output['Startx',i]+output['Endx',i])/2+tol, y=(output['Starty',i]+output['Endy',i])/2+tol, labels = output['T',i], col='grey', cex=1.5)
    }else if (slope<0){
      #
      # 兩點斜率為負時
      # --o--|   --o--|   --o--|
      #      o        |        |
      # -----|   -----o   ---o-|
      #
      # 不管甚麼狀況都要先畫從起點到(a,b)座標的那條平行線
      segments(output['Startx',i], output['Starty',i], a, b, col=output['Y',i], lwd = 10)
      if (output['Endx',i] == a & output['Endy',i] < b){
        # 補通過(a,b)、往(a,0)方向、長度<b的垂直線
        segments(a, b, a, output['Endy',i], col=output['Y',i], lwd = 10)
      }else if (output['Endx',i] == a & output['Endy',i] == 0){
        # 補通過(a,b)與(a,0)的垂直線
        segments(a, b, a, 0, col=output['Y',i], lwd = 10)
      }else if (output['Endx',i] < a & output['Endy',i] == 0){
        # 補通過(a,b)與(a,0)的垂直線
        segments(a, b, a, 0, col=output['Y',i], lwd = 10)
        # 補通過(a,0)、往(0,0)方向、長度<a的水平線
        segments(a, 0, output['Endx',i], 0, col=output['Y',i], lwd = 10)
      }
      text(x=(output['Startx',i]+output['Endx',i])/2, y=(output['Starty',i]+output['Endy',i])/2, labels = output['T',i], col='grey', cex=1.5)
    }else{
      #
      # 兩點斜率為正時
      # -----|   --o--|
      #      o        |
      # ---o-|   o----|
      #
      if (output['Startx',i] == a){
        segments(a, output['Starty',i], a, 0, col=output['Y',i], lwd = 10)
      }else if (output['Startx',i] < a){
        segments(output['Startx',i], output['Starty',i], a, b, col=output['Y',i], lwd = 10)
        segments(a, b, a, 0, col=output['Y',i], lwd = 10)
      }
      # 不管甚麼狀況最後都要畫到從(a,0)到終點座標的那條平行線
      segments(a, 0, output['Endx',i], 0, col=output['Y',i], lwd = 10)
      text(x=(output['Startx',i]+output['Endx',i])/2, y=(output['Starty',i]+output['Endy',i])/2, labels = output['T',i], col='grey', cex=1.5)
    }
  }
  points(c(output['Startx',],output['Endx',]), c(output['Starty',], output['Endy',]), type = 'p', lwd = 15, col='yellow')
}

#函數：走路流程的繪圖
addarrows = function(output, a, b){
  # 補arrows
  col_start_end_mend_all = {}
  for (i in 1:Wmax){
    
    output_temp = as.matrix(output[, which(output["Y",]==i)])
    output_order = order(output_temp["X",])
    col_start_end = {}
    
    for (j in 1:ncol(output_temp)){
      start = unname(output_temp[c("Startx","Starty"),output_order[j]])
      end = unname(output_temp[c("Endx","Endy"),output_order[j]])
      # 第一個位置，col_start_end第三欄1代表與上一筆有重複、反之則為0
      if (j == 1){
        col_start_end = unname(rbind(col_start_end, c(start,0)))
        col_start_end = unname(rbind(col_start_end, c(end,1)))
        # 第二個位置之後、且最後一個座標與新的座標重複時
      }else if (j > 1 & sum(col_start_end[nrow(col_start_end),1:2]==start)==2){
        col_start_end = unname(rbind(col_start_end, c(end,1)))
        # 第二個位置之後、且最後一個座標與新的座標無重複時
      }else if (j > 1 & sum(col_start_end[nrow(col_start_end),1:2]==start)!=2){
        col_start_end = unname(rbind(col_start_end, c(start,0)))
        col_start_end = unname(rbind(col_start_end, c(end,1)))
      }
    }
    # 補充U型生產線的轉折點(a,0)或(a,b)
    col_start_end_mend = {}
    for (k in 2:nrow(col_start_end)){
      if (col_start_end[k-1,1] == a & col_start_end[k,2] == 0 & col_start_end[k,3] == 1){
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k-1,])
        col_start_end_mend = rbind(col_start_end_mend, c(a,0,1))
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k,])
      } else if (col_start_end[k-1,2] == b & col_start_end[k,1] == a & col_start_end[k,3] == 1){
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k-1,])
        col_start_end_mend = rbind(col_start_end_mend, c(a,b,1))
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k,])
      } else {
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k-1,])
        col_start_end_mend = rbind(col_start_end_mend, col_start_end[k,])
      }
    }
    #移除重複
    col_start_end_mend = unique(col_start_end_mend)
    #補原點
    col_start_end_mend = rbind(col_start_end_mend, c(unname(output_temp[c("Startx","Starty"),output_order[1]]),0))
    #擴增欄位
    col_start_end_mend_aug = cbind(col_start_end_mend,0,0,0)
    for (k in 2:nrow(col_start_end_mend_aug)){
      #若不在U型路線上，計算步行距離/VW=步行時間,存於第五欄
      if (col_start_end_mend_aug[k,3]==0){
        col_start_end_mend_aug[k,5]=sqrt(sum((col_start_end_mend_aug[k-1,1:2]-col_start_end_mend_aug[k,1:2])^2))/VW
        #若在U型路線上，計算工作距離/VC=工作時間,存於第四欄
      }else{
        col_start_end_mend_aug[k,4]=sqrt(sum((col_start_end_mend_aug[k-1,1:2]-col_start_end_mend_aug[k,1:2])^2))/VC
      }
      #兩者加總存於第五欄
      col_start_end_mend_aug[k,6]=col_start_end_mend_aug[k,4]+col_start_end_mend_aug[k,5]
    }
    #補隸屬工作站
    col_start_end_mend_all = rbind(col_start_end_mend_all, cbind(col_start_end_mend_aug, rep(i, nrow(col_start_end_mend_aug))))
  }
  colnames(col_start_end_mend_all) = c('x','y','checked','working time','walking time','sum of time','workstation')
  col_start_end_mend_all
  
  for (l1 in 1:Wmax){
    temp = filter(as.data.frame(col_start_end_mend_all), workstation==l1)
    for (l2 in 2:nrow(temp)){
      
      if (temp$checked[l2]==0 & l2 == nrow(temp)){
        Arrows(temp$x[l2-1], temp$y[l2-1], temp$x[l2], temp$y[l2], lty = 'dotted', code = 2, arr.type="triangle", col = 'grey', arr.col = l1)
        text(x=(temp$x[l2-1]+temp$x[l2])/2, y=(temp$y[l2-1]+temp$y[l2])/2, labels = specify_decimal(temp$`sum of time`[l2],2), col = 'purple', cex = 1)
      }else if (temp$checked[l2]==0 & l2 != nrow(temp)){
        Arrows(temp$x[l2-1], temp$y[l2-1], temp$x[l2], temp$y[l2], lty = 'dotted', code = 2, arr.type="triangle", col = 'grey', arr.col = l1)
        text(x=(temp$x[l2-1]+temp$x[l2])/2, y=(temp$y[l2-1]+temp$y[l2])/2, labels = specify_decimal(temp$`sum of time`[l2],2), col = 'purple', cex = 1)
      }
    }
  }
}

#函數：將鄰近工作站的工作合併
merge = function(output, a, b){
  #計算每個工作站的總工作時間佔比、剃除只有一個工作的工作站
  Y_prob = filter(mutate(summarise(group_by(as.data.frame(t(output)),Y), freq = n(), sumT = sum(T)), prob = sumT/sum(sumT)), freq > 1)
  #依工作時間佔比抽工作站
  sampled_Y = sample(x=Y_prob$Y, size=1, prob = Y_prob$prob)
  
  #找出各工作站的邊邊
  temp = mutate(arrange(as.data.frame(t(output)),X), YY= abs(lag(Y)-Y) | abs(lead(Y)-Y))
  temp$YY[is.na(temp$YY)]=FALSE
  #選出位於邊邊的工作
  temp_filter = filter(temp, Y==sampled_Y, YY==TRUE)$index
  if (length(temp_filter)==1){
    sampled_X = temp_filter
  }else{
    sampled_X = sample(filter(temp, Y==sampled_Y, YY==TRUE)$index, size=1)
  }
  
  
  #若最後一個工作被選到
  
  index2X = temp[which(temp['index']==sampled_X),'X']
  
  if (index2X == nrow(temp)){
    temp[index2X, 'Y']=temp[index2X-1, 'Y']
    #若最第一個工作被選到
  } else if (index2X == 1){
    temp[index2X, 'Y']=temp[index2X+1, 'Y']
  } else{
    #若中間的工作被選到
    if (temp[index2X-1, 'YY']==TRUE & temp[index2X+1, 'YY']!=TRUE){
      temp[index2X, 'Y']=temp[index2X-1, 'Y']
    } else if (temp[index2X-1, 'YY']!=TRUE & temp[index2X+1, 'YY']==TRUE){
      temp[index2X, 'Y']=temp[index2X+1, 'Y']
    } else if (temp[index2X-1, 'YY']==TRUE & temp[index2X+1, 'YY']==TRUE){
      if (temp[index2X-1, 'Y']==sampled_Y & temp[index2X+1, 'Y']!=sampled_Y)
        temp[index2X, 'Y']=temp[index2X+1, 'Y']
      else if (temp[index2X-1, 'Y']!=sampled_Y & temp[index2X+1, 'Y']==sampled_Y){
        temp[index2X, 'Y']=temp[index2X-1, 'Y']
      }else{
        temp[index2X, 'Y']=sample(c(temp[index2X-1, 'Y'], temp[index2X+1, 'Y']), size=1)
      }
    }
  }
  
  output_change = t(arrange(as.data.frame(temp), index))
  
  return(output_change)
}






# 模擬退火法
E = sum(T)
# 疊代次數
H = 1000
# 連續不好的解之容忍次數
F = 25
# 模擬退火參數1
alpha = 0.1
# 模擬退火參數2
beta = 0.1
# 模擬退火參數3
Er = E
# 全域最佳解
Zglobal = {}
Zglobal$performance = E
# 連續不好的解的累計次數
f = 0
# 找搜尋解的總累計次數
h = 0
start.time <- Sys.time()
while (h < H){
  
  tryCatch({
    
    #要執行的指令放這裡
    O = solGen(Wmax, PR, T, a, b)
    while (length(unique(O['Y',])) != Wmax){
      O = solGen(Wmax, PR, T, a, b)
    }
    Z = evaluate(O, Wmax, a, b, VW, T)
    Zlocal = Z
    Onew = merge(O, a, b)
    Znew = evaluate(Onew, Wmax, a, b, VW, T)
    h = h + 1
    if (Znew$performance < Z$performance){
      O = Onew
      if (Znew$performance < Zlocal$performance){
        Olocal = O
        Zlocal = Znew
        f = 0
        f = f + 1
        if (f < F){
          E = E/(1+beta*E)
        }else{
          E = alpha*E + (1-alpha)*Er
          Er = E
          f = 0
        }
        if (Znew$performance < Zglobal$performance){
          Oglobal = O
          Zglobal = Znew
          f = 0
          f = f + 1
          if (f < F){
            E = E/(1+beta*E)
          }else{
            E = alpha*E + (1-alpha)*Er
            Er = E
            f = 0
          }
        }else{
          f = 0
          f = f + 1
          if (f < F){
            E = E/(1+beta*E)
          }else{
            E = alpha*E + (1-alpha)*Er
            Er = E
            f = 0
          }
        }
      }else{
        f = f + 1
        if (f < F){
          E = E/(1+beta*E)
        }else{
          E = alpha*E + (1-alpha)*Er
          Er = E
          f = 0
        }
      }
    }else{
      p = runif(n = 1)
      if (p < exp(-(Znew$performance-Z$performance)/E)){
        O = Onew
        f = f + 1
        if (f < F){
          E = E/(1+beta*E)
        }else{
          E = alpha*E + (1-alpha)*Er
          Er = E
          f = 0
        }
      }else{
        f = f + 1
        if (f < F){
          E = E/(1+beta*E)
        }else{
          E = alpha*E + (1-alpha)*Er
          Er = E
          f = 0
        }
      }
    }
    
  },warning = function(war){
    print(paste("MY_WARNING:  ",war)) #如果有warning則輸出warning,"MY_WARNING:  "這一行可以自己改
  },error = function(err) {
    print(paste("MY_ERROR:  ",err))   #如果有error則輸出error,"MY_Error:  "這一行可以自己改
  },finally = {
    print(paste(" ")) #最後一定要執行的指令或輸出
  })
  
  
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

plotU(Oglobal, a, b)
addarrows(Oglobal, a, b)
print(Oglobal[c('X','Y'),])
print(Zglobal)
print('1:黑；2：紅；3綠；4藍')

