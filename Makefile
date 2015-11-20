CC = gcc
CXX = g++

HOME=/home/yuguangyang/
ARMA_INCLUDE=-I$(HOME)Downloads/armadillo-5.100.2/include
RL_INCLUDE=-I../../include
GTEST_INCLUDE=-I$(HOME)workspace/libs/gtest-1.7.0/include
BOOST_INCLUDE=-I/opt/boost/boost_1_57_0
PROTO_INCLUDE=-I/usr/local/include

GTEST_PATH=-L$(HOME)workspace/libs/gtest-1.7.0/mybuilds
RL_PATH=-L../../lib
PROTO_PATH=-L/usr/local/lib

DEBUGFLAG=-DDEBUG -g3
RELEASEFLAG= -O3 -march=native -DARMA_NO_DEBUG
CXXFLAGS=  -std=c++0x -I$(HOME)Copy/workspace/munkres-cpp/src -D__LINUX  
#CXXFLAGS += $(DEBUGFLAG)
#CXXFLAGS += $(RELEASEFLAG)
LINKOPTFLAGS= -O3 -flto=4 -fwhole-program
LDFLAG=

OBJ=HungarianAlg.o controller.o model.o simulator.o Driver.o

test.exe: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAG) 
	

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(DEBUGFLAG) $^

clean:
	rm *.o *.exe
