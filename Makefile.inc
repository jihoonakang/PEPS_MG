# Main compiler with MPI. Use mpiifort for intel MPI.
FC = mpiifort

# Complie flag

# For debuggi#ng
FLAGS = -O0 -g -fPIC -check all -warn all -WB

# For moderate debugging
FLAGS = -O2 -fPIC -fp-model=precise -check bounds -fopenmp

# For optimized run
#FLAGS = -O3 -fPIC -fp-model=precise -xCORE-AVX512 -fopenmp
