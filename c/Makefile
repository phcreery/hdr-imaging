CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c99 -I.
LDFLAGS = -lm
USE_NML ?= 0

SRCS = main.c load_images.c gsolve.c hdr_debevec.c compute_irradiance.c tonemap.c

ifeq ($(USE_NML),1)
	CFLAGS += -DUSE_NML_LS -I./neat-matrix-library
	SRCS += neat-matrix-library/nml.c neat-matrix-library/nml_util.c
endif

OBJS = $(SRCS:.c=.o)

# Detect OS for output name
ifeq ($(OS),Windows_NT)
    TARGET = hdr.exe
else
    TARGET = hdr
endif

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJS) $(TARGET)

# Convenience target: build and run the test images
test: $(TARGET)
	./$(TARGET) ../images/1_4.png ../images/32_1.png output.png
