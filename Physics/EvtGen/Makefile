CXX = g++
CXXFLAGS = -g -Wall -fPIC 
LD = g++
LDFLAGS = -g -O2 -Wall -pthread -fPIC -shared 

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=libEvtGen.so

.SUFFIXES: .cpp

# ================================================================================
all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) -o libEvtGen.so

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $(SOURCES)
	
clean:
	rm -f libEvtGen.so  *.o 
