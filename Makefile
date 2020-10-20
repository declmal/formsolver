.PHONY: test clean

test:
	@mkdir -p build
	@cd build && cmake .. && make
	@mkdir -p bin
	@mv build/test_fem bin

clean:
	@rm -rf build bin
