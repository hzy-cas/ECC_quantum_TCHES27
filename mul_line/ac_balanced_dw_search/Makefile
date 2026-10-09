CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra
PROJECT_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CPPFLAGS += -DAC_BALANCED_SOURCE_DIR=\"$(PROJECT_ROOT)\"

BUILD_DIR := build
TARGET := $(BUILD_DIR)/ac_balanced_estimator
TEST_TARGET := $(BUILD_DIR)/test_schedule

CORE_SOURCES := \
	src/Infrastructure/GateManager.cpp \
	src/Infrastructure/MatrixLoader.cpp \
	src/BasicLogic/QuantumBasic.cpp \
	src/FieldArithmetic/QuantumField.cpp \
	src/FieldArithmetic/QuantumInversion.cpp \
	src/AlgorithmFlow/BalancedQuantumShor.cpp \
	src/AlgorithmFlow/BalancedResource.cpp

HEADERS := $(wildcard src/*/*.h) $(wildcard ../../circuit_backend/cpp/*)

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(CORE_SOURCES) src/main.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ src/main.cpp $(CORE_SOURCES)

$(TEST_TARGET): src/AlgorithmFlow/BalancedResource.cpp tests/test_schedule.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ tests/test_schedule.cpp src/AlgorithmFlow/BalancedResource.cpp

test: $(TARGET) $(TEST_TARGET)
	$(TEST_TARGET)
	python3 -m unittest discover -s tests -p 'test_*.py'

clean:
	rm -f $(TARGET) $(TEST_TARGET)

.PHONY: all test clean
