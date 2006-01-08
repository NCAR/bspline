
CCC=g++ -g -O0

all: lopass

OBJS = lopass.o flopass.o

FLIBS = -L/opt/SUNWspro/SC4.0/lib -lF77 -lM77 -lV77 -lsunmath
lopass : $(OBJS)
	$(CCC) -o lopass $(OBJS) $(FLIBS)


