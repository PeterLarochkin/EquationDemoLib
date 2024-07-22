n=125
threads_num=4
all: julia_experiment c_experiment c_omp_experiment julia_threads_experiment final
julia_experiment:
	julia --project=. ./src/Main.jl $n $n > ./docs/results/julia.txt
julia_threads_experiment:
	export JULIA_NUM_THREADS=${threads_num} && julia -t ${threads_num} --project=. ./src/MainParallel.jl $n $n > ./docs/results/julia_parallel.txt
c_experiment:
	(gcc-14 -o ./c_analogue/proximity_search ./c_analogue/proximity_search.c && ./c_analogue/proximity_search $n $n) > ./docs/results/c.txt && rm ./c_analogue/proximity_search
c_omp_experiment:
	(gcc-14 -fopenmp ./c_analogue/proximity_search_OMP.c -o ./c_analogue/proximity_search_OMP && export OMP_NUM_THREADS=${threads_num} && ./c_analogue/proximity_search_OMP $n $n) > ./docs/results/c_OMP.txt && rm ./c_analogue/proximity_search_OMP
generate_docs:
	cd ./docs && julia -e 'using Pkg; Pkg.instantiate()' && julia -e 'using Pkg; Pkg.add("Documenter")' && julia make.jl && cd - && echo "see file ./docs/build/index.html"
final:
	echo "\nJulia single:" && tail -n 1 ./docs/results/julia.txt && \
	echo "\nJulia parallel:" && tail -n 1 ./docs/results/julia_parallel.txt && \
	echo "\nC:" && tail -n 1 ./docs/results/c.txt && \
	echo "\nC openMP :" && tail -n 1 ./docs/results/c_OMP.txt
