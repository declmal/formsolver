.PHONY: form ccx test clean

ifeq ($(OS),Windows_NT)
  SHARED_LIBRARY_SUFFIX := dll
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S), Darwin)
    SHARED_LIBRARY_SUFFIX := dylib
  else
    SHARED_LIBRARY_SUFFIX := so
  endif
endif

form:
	@mkdir -p build
	@cd build && cmake .. && make
	@mkdir -p lib
	@mv build/libform.$(SHARED_LIBRARY_SUFFIX) lib

ccx:
	@mkdir -p build/ccx
	@cd build/ccx && cmake ../../ccx && make
	@mkdir -p lib
	@mv build/ccx/libccx.${SHARED_LIBRARY_SUFFIX} lib

test: form ccx
	@mkdir -p build/tests
	@cd build/tests && cmake ../../tests && make
	@mkdir -p bin
	@mv build/tests/test_* bin

clean:
	@rm -rf build lib bin
