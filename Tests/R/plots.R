
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
for (cutoff in c(5, 30, 0))
{
    for (bc in c(0, 1, 2))
    {

        output = paste("spline-", cutoff, "-", bc, ".out", sep="")
        cmd = paste("../C++/bspline -w", cutoff, "-b", bc, "-i", filename, "-o", output)
        cat(cmd, "\n")
        system(cmd)
        spline = read.table(output, header=TRUE)
        plot(input[,1], input[,2], type="p", main="Input Temperatures",
            xlab="Seconds", ylab="Temp (C)")
        plot(spline[,1], spline[,3], type="l",
            main=paste("Smoothed Temperatures: bc=",bc,"cutoff=",cutoff),
            xlab="Seconds", ylab="Temp (C)")
    }
}

graphics.off()
