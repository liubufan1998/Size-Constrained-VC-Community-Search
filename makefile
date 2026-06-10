# 编译器和参数
CXX := g++
CXXFLAGS := -std=c++11 -O3 -flto

# 目标文件
TARGET := SVC-TD
SRCS := main.cpp Graph.cpp VC_tool.cpp
OBJS := $(SRCS:.cpp=.o) 

# 默认目标
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

clean-all:
	find . -type f \( -name 'CMakeFiles' -o -name 'cmake_install.cmake' -o -name 'CMakeCache.txt' -o -name 'Makefile' -o -name '*.cmake' \) -delete
