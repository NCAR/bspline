

generate_plot <- function(fieldname, input, cutoff, bc)
{
    # remove initial offset from times
    input[,1] <- input[,1] - input[1,1]
    basename = paste("spline-", fieldname, "-", cutoff, "-", bc, sep="")
    output = paste(basename, ".out", sep="")
    filename = paste("../Data/sample.", fieldname, sep="")
    cmd = paste("../C++/bspline -w", cutoff, "-b", bc, "-i", filename, "-o", output)
    cat(cmd, "\n")
    system(cmd)
    spline = read.table(output, header=TRUE)
    png(paste("plot-", basename, ".png", sep=""), width=1000, height=800)
    par(mfrow=c(2,1))
    plot(input[,1], input[,2], type="p", main="Input Temperatures",
         xlab="Seconds", ylab="Temp (C)")
    # png(paste("plot-", basename, ".png"), width=1000, height=800)
    plot(spline[,1], spline[,3], type="l",
         main=paste("Smoothed Temperatures: bc=",bc,"cutoff=",cutoff),
         xlab="Seconds", ylab="Temp (C)")
    graphics.off()
}

temps = matrix(scan("../Data/sample.temps"), ncol=2, byrow=1, 
               dimnames=list(c(), c("time","temp")))
wdir =  matrix(scan("../Data/sample.wdir"), ncol=2, byrow=1, 
               dimnames=list(c(), c("time","wdir")))
wspd =  matrix(scan("../Data/sample.wspd"), ncol=2, byrow=1, 
               dimnames=list(c(), c("time","wspd")))


graphics.off()
# postscript("plots.ps")
# x11()
for (cutoff in c(5, 30, 0))
{
    for (bc in c(0, 1, 2))
    {
        generate_plot("temps", temps, cutoff, bc)
        generate_plot("wdir", wdir, cutoff, bc)
        generate_plot("wspd", wspd, cutoff, bc)
    }
}

