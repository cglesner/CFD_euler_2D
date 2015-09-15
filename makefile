Run_Simulation: all
	g++ -Wall -o Run_Simulation main.o Euler_Helper_Functions.o        \
	Roe_Flux.o Van_Leer_Flux.o Geometry.o Time_Step.o simulation_io.o  \
	diagnostics.o mms.o muscl.o boundary_conditions.o -g

all: main.o Euler_Helper_Functions.o Roe_Flux.o Van_Leer_Flux.o Geometry.o \
	Time_Step.o simulation_io.o mms.o muscl.o boundary_conditions.o diagnostics.o

main.o: main.cpp main.h Geometry.h Time_Step.h Roe_Flux.h Van_Leer_Flux.h  \
	simulation_io.h Euler_Helper_Functions.h mms.h muscl.h boundary_conditions.h
	g++ -Wall -c main.cpp -g

Euler_Helper_Functions.o: Euler_Helper_Functions.cpp Euler_Helper_Functions.h
	g++ -Wall -c Euler_Helper_Functions.cpp -g

Roe_Flux.o: Roe_Flux.cpp Roe_Flux.h Euler_Helper_Functions.h
	g++ -Wall -c Roe_Flux.cpp -g

Van_Leer_Flux.o: Van_Leer_Flux.cpp Van_Leer_Flux.h Euler_Helper_Functions.h
	g++ -Wall -c Van_Leer_Flux.cpp -g

Geometry.o: Geometry.cpp Geometry.h
	g++ -Wall -c Geometry.cpp -g

Time_Step.o: Time_Step.cpp Time_Step.h
	g++ -Wall -c Time_Step.cpp -g

Flux_Determiner.o: Flux_Determiner.cpp Flux_Determiner.h
	g++ -Wall -c Flux_Determiner.cpp -g

Roe_Flux_Determiner.o: Roe_Flux_Determiner.cpp Roe_Flux_Determiner.h
	g++ -Wall -c Roe_Flux_Determiner.cpp -g

simulation_io.o: simulation_io.cpp simulation_io.h
	g++ -Wall -c simulation_io.cpp -g

mms.o: mms.cpp mms.h
	g++ -Wall -c mms.cpp -g

muscl.o: muscl.cpp muscl.h
	g++ -Wall -c muscl.cpp -g

boundary_conditions.o: boundary_conditions.cpp boundary_conditions.h
	g++ -Wall -c boundary_conditions.cpp -g

diagnostics.o: diagnostics.cpp diagnostics.h
	g++ -Wall -c diagnostics.cpp -g

#TestExec: Unit_Tests.cpp all
#	g++ -Wall -o TestExec Unit_Tests.cpp Geometry.o Time_Step.o \
#	Van_Leer_Flux.o Euler_Helper_Functions.o Roe_Flux.o         \
#	Flux_Determiner.o Roe_Flux_Determiner.o simulation_io.o mms.o muscl.o -g

clean:
	rm *.o TestExec Run_Simulation



