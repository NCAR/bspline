
temps = matrix(scan("sample.temps"), ncol=2, byrow=1, 
	dimnames=list(c(), c("time","temp")))

wdir =  matrix(scan("sample.wdir"), ncol=2, byrow=1, 
	dimnames=list(c(), c("time","wdir")))

wspd =  matrix(scan("sample.wspd"), ncol=2, byrow=1, 
	dimnames=list(c(), c("time","wspd")))

filename = "sample.temps"
temps[,1] <- temps[,1] - temps[1,1]
input = temps 

graphics.off()
# postscript("plots.ps")
png("plot-%d.png", width=600, height=800)
# x11()
par(mfrow=c(2,1))
for (cutoff in c(5, 30, 0)) {
for (bc in c(0, 1, 2)) {

	cmd = paste("../C++/driver 1", cutoff, bc, "<", filename)
	cat(cmd, "\n")
	system(cmd)
	output = paste("spline-", cutoff, "-", bc, ".out", sep="")
	system(paste("mv -f spline.out", output))
	spline = matrix(scan(output), ncol=3, byrow=1,
		dimnames=list(c(), c("time","smoothed","deriv")))
	input = matrix(scan("input.out"), ncol=3, byrow=1,
		dimnames=list(c(), c("time","input","smoothed")))
	plot(input[,1], input[,2], type="p", main="Input Temperatures",
		xlab="Seconds", ylab="Temp (C)")
	plot(spline[,1], spline[,2], type="l",
		main=paste("Smoothed Temperatures: bc=",bc,"cutoff=",cutoff),
		xlab="Seconds", ylab="Temp (C)")
}
}

graphics.off()
