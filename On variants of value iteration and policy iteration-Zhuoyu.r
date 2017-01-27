# On variants of value iteration and policy iteration:

# Modified policy iteration: instead of fully (to convergence) evaluating a policy, just compute a few Bellman backups and then improve. See section 6.5 of Puterman (1994)
# Asynchronous value iteration, Gauss-Seidel and Jacobi variants. See section 6.3.3 of Puterman (1994)
# Andrew W. Moore and Christopher G. Atkeson. Prioritized sweeping: Reinforcement learning with less data and less real time. Machine Learning, 13, 1993.


# Example 4.1 on page 82 in Sutton 2016.

# demo of the policy evaluation on page 83 of (Sutton 2016)
# iterative policy evaluation - two array version
# also called Jacobi-style algorithm
# 6 instead of 4 is used here: the top, bottom, left and right lines/columns are temp values, the 4*4 in the middle are the matrix.
V0 <- matrix(0, nrow = 6, ncol = 6)  # V(s)
# random policy
policy <- array(0.25, dim=c(4,4,4)) # in order of up, right, down, left,

gamma <- 1 # discount rate
r <- -1 	# reward 
delta=100
k=-1
theta <- 0.001 # stopping criteria
while(delta>=theta){
	k <- k + 1
	delta<- 0
	V1 <- V0
	cat("k=",k,"\n")
	print(V1[2:5, 2:5])
	cat("\n")
	
	V0[2:5, 2:5] <- policy[,,1]*(r+gamma*V1[1:4, 2:5]) +
			policy[,,2]*(r+gamma*V1[2:5, 3:6]) +
			policy[,,3]*(r+gamma*V1[3:6, 2:5]) +
			policy[,,4]*(r+gamma*V1[2:5, 1:4])
	V0[2,2] <- 0
	V0[5,5] <- 0
	V0[1, 2:5] <- V0[2, 2:5]
	V0[6, 2:5] <- V0[5, 2:5]
	V0[2:5, 1] <- V0[2:5, 2]
	V0[2:5, 6] <- V0[2:5, 5]
	delta <- max(delta, abs(V1[2:5, 2:5]-V0[2:5, 2:5]))
}
cat("Reach Convergence at k=", k, ", with the value function of", "\n")
	print(V1[2:5, 2:5])
	cat("\n")
	

#----------------------------------------------------------------------------------------
# iterative policy evaluation - in place version
# Also called as Gauss-Seidel-style algorithm
# 6 instead of 4 is used here: the top, bottom, left and right lines/columns are temp values, the 4*4 in the middle are the matrix.
V0 <- matrix(0, nrow = 6, ncol = 6)  # V(s)
# random policy
policy <- array(0.25, dim=c(4,4,4)) # in order of up, right, down, left,

gamma <- 1 # discount rate
r <- -1 	# reward 
delta=100
k=-1
theta <- 0.001 # stopping criteria
while(delta>=theta){
	k <- k + 1
	delta<- 0
	V1 <- V0
	cat("k=",k,"\n")
	print(V1[2:5, 2:5])
	cat("\n")
	
	for (i in 2:5){
		for(j in 2:5){
			if(!((i==2&j==2)|(i==5&j==5))){
				V0[i,j] <- policy[i-1,j-1,1]*(r+gamma*V0[i-1, j]) +
						policy[i-1,j-1,2]*(r+gamma*V0[i, j+1]) +
						policy[i-1,j-1,3]*(r+gamma*V0[i+1, j]) +
						policy[i-1,j-1,4]*(r+gamma*V0[i, j-1])
			}
		}
	}
	V0[1, 2:5] <- V0[2, 2:5]
	V0[6, 2:5] <- V0[5, 2:5]
	V0[2:5, 1] <- V0[2:5, 2]
	V0[2:5, 6] <- V0[2:5, 5]
	delta <- max(delta, abs(V1[2:5, 2:5]-V0[2:5, 2:5]))
}
cat("Reach Convergence at k=", k, ", with the value function of", "\n")
	print(V1[2:5, 2:5])
	cat("\n")



#-----------------------------------------------------------------------------------------------------------
# Prioritized Sweeping By Moore et al., 1993
# 6 instead of 4 is used here: the top, bottom, left and right lines/columns are temp values, the 4*4 in the middle are the matrix.
V0 <- matrix(0, nrow = 6, ncol = 6)  # V(s)

# random policy
policy <- array(0.25, dim=c(4,4,4)) # in order of up, right, down, left,

gamma <- 1 # discount rate
r <- -1 	# reward 
delta=100
k=-1
theta <- 0.001 # stopping criteria
eplison <- 0.001 # priority criteria
priority_queue <- 2
priority_p <- 0.1
while((length(priority_queue)>0)|(delta>=theta)){
	k <- k + 1
	delta<- 0
	V1 <- V0
	cat("k=",k,"\n")
	print(V1[2:5, 2:5])
	cat("\n")
	i <- priority_queue[1]%/%4 + 2
	j <- priority_queue[1]%%4 + 1
	if(j==1){
		j <- 5
		i <- i-1
	}
	if(length(priority_queue)>1){
		priority_queue <- priority_queue[2:length(priority_queue)]
		priority_p <- priority_p[2:length(priority_p)]
	} else {
		priority_queue <- NULL
		priority_p <- NULL
	}
	V0[i, j] <- policy[i-1, j-1,1]*(r+gamma*V1[i-1, j]) +
			policy[i-1, j-1,2]*(r+gamma*V1[i, j+1]) +
			policy[i-1, j-1,3]*(r+gamma*V1[i+1, j]) +
			policy[i-1, j-1,4]*(r+gamma*V1[i, j-1])

	new_delta <- abs(V1[i,j]-V0[i,j])
	# V1[i,j] <- V0[i,j]
	if(i==2){
		V0[1,j] <- V0[i,j]
	}
	if(i==5){
		V0[6,j] <- V0[i,j]	
	}
	if(j==2){
		V0[i,1] <- V0[i,j]	
	}
	if(j==5){
		V0[i,6] <- V0[i,j]	
	}
	V0[2,2] <-0
	V0[5,5] <-0	
	delta <- max(delta, new_delta)
	
	if(new_delta>eplison){
		x <- c(max(i-1,2), i, min(i+1,5), i)
		y <- c(j, min(j+1,5), j, max(j-1,2))
		index <- (x-2)*4+(y-1)
		for(m in index)
		{
			if(m!=1& m!=16){
				if(m %in% priority_queue){
					if(priority_p[which(priority_queue==m)]<new_delta){
						priority_p[which(priority_queue==m)] <- new_delta
						new_order <- order(priority_p,decreasing=TRUE)
						priority_p <- priority_p[new_order]
						priority_queue <- priority_queue[new_order]
					}
				} else {
					priority_queue <- c(priority_queue,m)
					priority_p <- c(priority_p, new_delta)
					new_order <- order(priority_p,decreasing=TRUE)
					priority_p <- priority_p[new_order]
					priority_queue <- priority_queue[new_order]
				}
			}
		}
	}
	
}
cat("Reach Convergence at k=", k, ", with the value function of", "\n")
	print(V1[2:5, 2:5])
	cat("\n")

	
#-------------------------------------------------------------------------------------------------		
#  Policy iteration (using iterative policy evaluation) on page 87 of (Sutton, 2016)
# 6 instead of 4 is used here: the top, bottom, left and right lines/columns are temp values, the 4*4 in the middle are the matrix.
V0 <- matrix(0, nrow = 6, ncol = 6)  # V(s)
# initial policy
# "u","r","d","l": in order of up, right, down, left,
policy <- array(0.25, dim=c(4,4,4)) # in order of up, right, down, left,
policy[1,1,] <-0
policy[4,4,] <-0
#policy0 <- matrix(sample(c("u","r","d","l"), 4*4, replace=TRUE), nrow=4, ncol=4)
# policy0[1,1] <- "ter"
# policy0[4,4] <- "ter"
gamma <- 1 # discount rate
r <- -1 	# reward 

delta=100
delta2=100
k1=-1
k2=-1
theta <- 0.00001 # stopping criteria for 
theta2 <- 0.00001 # 
while(delta2>=theta2) {
	k2 <- k2 + 1
	# policy evaluation step
	while(delta>=theta){
		k1 <- k1 + 1
		delta<- 0
		V1 <- V0
	#	cat("k=",k,"\n")
	#	print(V1[2:5, 2:5])
	#	cat("\n")
		
		for (i in 2:5){
			for(j in 2:5){
				if(!((i==2&j==2)|(i==5&j==5))){
					V0[i,j] <- policy[i-1,j-1,1]*(r+gamma*V0[i-1, j]) +
								policy[i-1,j-1,2]*(r+gamma*V0[i, j+1]) +
								policy[i-1,j-1,3]*(r+gamma*V0[i+1, j]) +
								policy[i-1,j-1,4]*(r+gamma*V0[i, j-1])
				}
			}
		}
		V0[1, 2:5] <- V0[2, 2:5]
		V0[6, 2:5] <- V0[5, 2:5]
		V0[2:5, 1] <- V0[2:5, 2]
		V0[2:5, 6] <- V0[2:5, 5]
		delta <- max(delta, abs(V1[2:5, 2:5]-V0[2:5, 2:5]))
	}
	# cat("Reach Convergence at k1=", k1, ", with the value function of", "\n")
	#	print(V1[2:5, 2:5])
	#	cat("\n")

	# Policy improvement step
	# policy_stable <- TRUE
	old_actions <- policy
	for (i in 1:4){
		for(j in 1:4){
			if(!((i==1&j==1)|(i==4&j==4))){
				temp_v <- c(V0[i, j+1], V0[i+1, j+2], V0[i+2, j+1], V0[i+1, j])
				argmax <- which(temp_v==max(temp_v))
				policy[i,j,] <-0
				policy[i,j,argmax] <- 1/length(argmax)
			}
		}
	}
	delta2<-max(abs(old_actions-policy))
	print(policy)
}
print(V1[2:5, 2:5])
cat("\n")
print(policy)
print(k1)
print(k2)




#-----------------------------------------------------------------------------------------
# value iteration . Section 4.4 of Sutton 2016 (page 89)
# It is also called as modified policy iteration in Puterman (1994)
# 6 instead of 4 is used here: the top, bottom, left and right lines/columns are temp values, the 4*4 in the middle are the matrix.
V0 <- matrix(0, nrow = 6, ncol = 6)  # V(s)
# initial policy
# "u","r","d","l": in order of up, right, down, left,
policy <- array(0.25, dim=c(4,4,4)) # in order of up, right, down, left,
policy[1,1,] <-0
policy[4,4,] <-0
#policy0 <- matrix(sample(c("u","r","d","l"), 4*4, replace=TRUE), nrow=4, ncol=4)
# policy0[1,1] <- "ter"
# policy0[4,4] <- "ter"
gamma <- 1 # discount rate
r <- -1 	# reward 

delta=100

k1=-1

theta <- 0.00001 # stopping criteria for 

while(delta>=theta){
	k1 <- k1 + 1
	delta<- 0
	V1 <- V0
#	cat("k=",k,"\n")
	print(V1[2:5, 2:5])
#	cat("\n")
	
	for (i in 2:5){
		for(j in 2:5){
			if(!((i==2&j==2)|(i==5&j==5))){
				V0[i,j] <- max(r+gamma*V1[i-1, j],r+gamma*V1[i, j+1],r+gamma*V1[i+1, j],r+gamma*V1[i, j-1])
			}
		}
	}
	V0[1, 2:5] <- V0[2, 2:5]
	V0[6, 2:5] <- V0[5, 2:5]
	V0[2:5, 1] <- V0[2:5, 2]
	V0[2:5, 6] <- V0[2:5, 5]
	delta <- max(delta, abs(V1[2:5, 2:5]-V0[2:5, 2:5]))
}
# cat("Reach Convergence at k1=", k1, ", with the value function of", "\n")
print(V1[2:5, 2:5])
#	cat("\n")


for (i in 1:4){
	for(j in 1:4){
		if(!((i==1&j==1)|(i==4&j==4))){
			temp_v <- c(V0[i, j+1], V0[i+1, j+2], V0[i+2, j+1], V0[i+1, j])
			argmax <- which(temp_v==max(temp_v))
			policy[i,j,] <-0
			policy[i,j,argmax] <- 1/length(argmax)
		}
	}
}


print(V1[2:5, 2:5])
cat("\n")
print(policy)
print(k1)

