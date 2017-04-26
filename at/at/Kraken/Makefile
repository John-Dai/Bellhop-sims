
EXECUTABLES = bounce.exe field.exe field3d.exe kraken.exe krakenc.exe modasc.exe modbin.exe

BOUNCE_OBJECTS = bounce.o BCImpedancec.o KrakencMod.o kuping.o PekerisRoot.o ../misc/libmisc.a

BOUNCE_LIBS = -L../misc -lmisc

FIELD_OBJECTS = field.o evaluate.o evaluatead.o evaluatecm.o getmodes.o PekerisRoot.o ../misc/libmisc.a

FIELD_LIBS = -L../misc -lmisc

FIELD3D_OBJECTS = field3d.o ElementMod.o evaluate3d.o evaluategb.o evaluatepdq.o getmodes.o PekerisRoot.o ../misc/libmisc.a

FIELD3D_LIBS = -L../misc -lmisc

KRAKEN_OBJECTS = BCImpedance.o KrakenMod.o kuping.o PekerisRoot.o InverseIterationMod.o RootFinderSecantMod.o zbrentx.o kraken.o ../misc/libmisc.a

KRAKEN_LIBS = -L../misc -lmisc

KRAKENC_OBJECTS = BCImpedancec.o KrakencMod.o kuping.o PekerisRoot.o InverseIterationMod.o RootFinderSecantMod.o krakenc.o ../misc/libmisc.a

KRAKENC_LIBS = -L../misc -lmisc

MODASC_OBJECTS = modasc.o

MODASC_LIBS = 

MODBIN_OBJECTS = modbin.o

MODBIN_LIBS = 

all: $(EXECUTABLES)
	@echo "KRAKEN built"
	@echo "************"
	@echo " "

install: all
	for f in $(EXECUTABLES) ; do \
		echo "----- Installing $$f"; cp -p $$f ../bin; \
	done

clean:
	-rm -f *.o *.mod *.exe *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 

# Declare explicit dependencies (and explicit build rules as needed)

bounce.exe: $(BOUNCE_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BOUNCE_OBJECTS) $(BOUNCE_LIBS)

field.exe: $(FIELD_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELD_OBJECTS) $(FIELD_LIBS)

field3d.exe: $(FIELD3D_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELD3D_OBJECTS) $(FIELD3D_LIBS)

kraken.exe: $(KRAKEN_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(KRAKEN_OBJECTS) $(KRAKEN_LIBS)

krakenc.exe: $(KRAKENC_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(KRAKENC_OBJECTS) $(KRAKENC_LIBS)

modasc.exe: $(MODASC_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(MODASC_OBJECTS) $(MODASC_LIBS)

modbin.exe: $(MODBIN_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(MODBIN_OBJECTS) $(MODBIN_LIBS)

ElementMod.o elementmod.mod: ElementMod.f90
	$(FC) -c $(FFLAGS) $<

InverseIterationMod.o inverseiterationmod.mod: InverseIterationMod.f90
	$(FC) -c $(FFLAGS) $<

RootFinderSecantMod.o rootfindersecantmod.mod: RootFinderSecantMod.f90
	$(FC) -c $(FFLAGS) $< 

KrakenMod.o krakenmod.mod: KrakenMod.f90
	$(FC) -c $(FFLAGS) $< 

KrakencMod.o krakencmod.mod: KrakencMod.f90
	$(FC) -c $(FFLAGS) $< 

BCImpedance.o: krakenmod.mod BCImpedance.f90
BCImpedancec.o: krakencmod.mod BCImpedancec.f90
bounce.o: krakencmod.mod bounce.f90
covar.o: covar.f90
evaluate.o: evaluate.f90
evaluate3d.o: elementmod.mod evaluate3d.f90
evaluatead.o: evaluatead.f90
evaluatecm.o: evaluatecm.f90
evaluategb.o: elementmod.mod evaluategb.f90
evaluatepdq.o: elementmod.mod evaluatepdq.f90
field.o: field.f90
field3d.o: elementmod.mod field3d.f90
getmodes.o: getmodes.f90
kraken.o: krakenmod.mod kraken.f90
krakenc.o: krakencmod.mod krakenc.f90
kuping.o: kuping.f90
modasc.o: modasc.f90
modbin.o: modbin.f90
PekerisRoot.o: PekerisRoot.f90
zbrentx.o: zbrentx.f90
zsecantcx.o: zsecantcx.f90
zsecantx.o: zsecantx.f90

