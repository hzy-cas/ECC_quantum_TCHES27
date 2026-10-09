CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall
TARGET = shor_estimate

SRCDIR = src
SOURCES = $(SRCDIR)/main.cpp \
          $(SRCDIR)/Infrastructure/GateManager.cpp \
          $(SRCDIR)/Infrastructure/MatrixLoader.cpp \
          $(SRCDIR)/BasicLogic/QuantumBasic.cpp \
          $(SRCDIR)/FieldArithmetic/QuantumField.cpp \
          $(SRCDIR)/FieldArithmetic/QuantumInversion.cpp \
          $(SRCDIR)/AlgorithmFlow/QuantumShor.cpp

all: $(TARGET)

$(TARGET): $(SOURCES) $(wildcard ../../circuit_backend/cpp/*)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
