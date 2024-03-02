JL = julia --project=@

default: format

format:
	$(JL) -e 'using JuliaFormatter; format("src", BlueStyle())'
