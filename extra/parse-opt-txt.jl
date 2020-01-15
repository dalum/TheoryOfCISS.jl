function parse_opt_file(filename)
    lines = readlines(filename)
    idx1 = findlast(line -> occursin("Band", line) & occursin("Eigenvalue", line), lines)
    idx2 = findfirst(line -> occursin("Gap", line), lines[idx1+1:end])
    return map(line -> parse(Float64, getindex(split(line), 2)), lines[idx1+1:idx1+idx2-2])
end
