
TARGETS = benchmark_math benchmark_rrc mean_free_path

all: $(TARGETS)

benchmark_math: benchmark_math.c ../util/util_misc.c
	$(CC) -g -O2 -pthread -Wall -I../util -I.. -lrt -lm -o $@ $^

benchmark_rrc: benchmark_rrc.c ../util/util_misc.c
	$(CC) -g -O2 -pthread -Wall -I../util -I.. -lrt -lm -o $@ $^

mean_free_path: mean_free_path.c ../util/util_misc.c
	$(CC) -g -O2 -pthread -Wall -I../util -I.. -lrt -lm -o $@ $^

clean:
	rm -f $(TARGETS)

