---
title: 国庆画个国旗
tags: coding
---

使用R，参考[这篇文章](https://blog.csdn.net/qq_32832803/article/details/82908259)。

```R
star <- function(posXY, size=1, theta=0, color="yellow")
{
	alpha <- 2 * pi / 5
	R <- size
	r <- sin(pi/10) * R / sin(7*pi/10)
	pointpos <- matrix(0, nrow=10, ncol=2)
	Beta <- (pi/10+theta) + (0:4) * alpha
	pointpos[seq(1,9,2),] <- cbind(cos(Beta)*R+posXY[1], sin(Beta)*R+posXY[2])
	Beta <- (3*pi/10+theta) + (0:4) * alpha
	pointpos[seq(2,10,2),] <- cbind(cos(Beta)*r+posXY[1], sin(Beta)*r+posXY[2])
	polygon(pointpos, col=color, border=color)
}

rotTheta <- function(posxy, posXY)
{
	atan((posXY[2]-posxy[2])/(posXY[1]-posxy[1])) + pi / 2
}


width <- 3
height <- 2 * width / 3
d <- width / 30
yshift <- height / 2

Flag <- cbind(c(0, width, width, 0), c(0, 0, height, height))
posXY <- c(5*d, yshift+5*d)
posxy <- matrix(c(10, 8, 12, 6, 12, 3, 10, 1), ncol=2, byrow=T)*d+matrix(c(0, yshift), nrow=4, ncol=2, byrow=T)
rottheta <- rep(0, nrow(posxy))
rottheta <- apply(posxy, 1, rotTheta, posXY)

plot(Flag, type="n", col="red", axes=F, xlim=c(0,width), ylim=c(0, height), xlab="", ylab="", main="People's Republic of China", asp=1)
polygon(Flag, density=NULL, col="red", border="red")
star(posXY, size=3*height/20, theta=0, color="yellow")

for(i in 1:nrow(posxy))
{
	star(posxy[i,], size=height/20, theta=rottheta[i], color="yellow")
}
```


![loveCN](https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master/downloads/images/national-flag-cn.jpeg)

祝祖国母亲生日快乐。