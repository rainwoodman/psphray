FSRC = rwarn.f daxpy.f dscal.f d1mach.f idamax.f ddot.f bnorm.f cfode.f dgbfa.f dgbsl.f dgefa.f dgesl.f ewset.f fdump.f fnorm.f i1mach.f intdy.f j4save.f lsoda.f prja.f solsy.f stoda.f vmnorm.f

all: main.o liblsoda.a
	gcc -g -O0 -o a.out $^ -lm -lgfortran

liblsoda.a: $(FSRC:.f=.o) lsoda_wrapper.o
	ar cur $@ $^

#lsoda_combined.f: $(FSRC)
#	cat $^ > $@

.f.o:
	gfortran -g -O0 -c -o $@ $<

.c.o:
	gcc -g -O0 -c -o $@ $<


