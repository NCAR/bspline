
temps = matrix(scan("sample.temps"), ncol=2, byrow=1, 
	dimnames=list(c(), c("time","temp")))

wdir =  matrix(scan("sample.wdir"), ncol=2, byrow=1, 
	dimnames=list(c(), c("time","wdir")))

wspd =  matrix(scan("sample.wspd"), ncol=2, byrow=1, 
	dimnames=list(c(), c("time","wspd")))

filename = "sample.temps"
temps[,1] <- temps[,1] - temps[1,1]
input = temps 

postscript("plots.ps")
par(mfrow=c(7,1))
plot(input[,1], input[,2], type="p", main="Input Temperatures",
	xlab="Seconds", ylab="Temperature (Celsisus)")

for (cutoff in c(5, 30)) {
for (bc in c(0, 1, 2)) {

	cmd = paste("driver 1", cutoff, bc, "<", filename)
	cat(cmd, "\n")
	system(cmd)
	spline = matrix(scan("spline.out"), ncol=3, byrow=1,
		dimnames=list(c(), c("time","smoothed","deriv")))
	input = matrix(scan("input.out"), ncol=3, byrow=1,
		dimnames=list(c(), c("time","input","smoothed")))
	plot(spline[,1], spline[,2], type="p",
		main=paste("Smoothed Temperatures: bc=",bc,"cutoff=",cutoff),
		xlab="Seconds", ylab="Temperature (Celsisus)")
}
}

