################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AlignmentRecord.cpp \
../src/BandAligner.cpp \
../src/Config.cpp \
../src/HashTableLongKey.cpp \
../src/HashTableMethod.cpp \
../src/HashTableShortKey.cpp \
../src/QueryDataset.cpp \
../src/QueryRead.cpp \
../src/SubjectDataset.cpp \
../src/SubjectRead.cpp \
../src/SubjectReadAligner.cpp \
../src/main.cpp 

OBJS += \
./src/AlignmentRecord.o \
./src/BandAligner.o \
./src/Config.o \
./src/HashTableLongKey.o \
./src/HashTableMethod.o \
./src/HashTableShortKey.o \
./src/QueryDataset.o \
./src/QueryRead.o \
./src/SubjectDataset.o \
./src/SubjectRead.o \
./src/SubjectReadAligner.o \
./src/main.o 

CPP_DEPS += \
./src/AlignmentRecord.d \
./src/BandAligner.d \
./src/Config.d \
./src/HashTableLongKey.d \
./src/HashTableMethod.d \
./src/HashTableShortKey.d \
./src/QueryDataset.d \
./src/QueryRead.d \
./src/SubjectDataset.d \
./src/SubjectRead.d \
./src/SubjectReadAligner.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -fopenmp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


