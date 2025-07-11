module I_nr
using CSV, SpecialFunctions

export I

inr_table = Vector{Vector{Float32}}()  # Vector of vectors (matrix-like)

# Open the CSV file and process it line by line
open("../utils/inr_output_250_250.csv", "r") do file
	# Iterate over each line of the file
	for line in eachline(file)
		# Split the line by commas (CSV format) and push it as a new vector to the matrix
		row = [parse(Float32, d) for d in split(line, ",")[1:end-1]]
		push!(inr_table, row)
	end
end


find_row(n, r) = div(n*(n+1),2) + 1 + r # index of the vector containing values for I_{n, r}


"""
    nearest_grid_idx(n, r, x)

Find the stored value of ``I_{n, r}`` closest to the supplied argument ``x``.
"""
function nearest_grid_idx(n, r, x)
	# get logarithms of all the values, since the precomputed x-values are on a log scale
	log_xmin = log10(max(10^-4, 0.01 * (sqrt(n) - sqrt(r))^2))
    log_xmax = log10(10 * (sqrt(n+1) + sqrt(r+1))^2)
	log_x = log10(x)

	if log_x > log_xmax return -1; end
	if log_x < log_xmin return 0; end 
	
	return (log_x-log_xmin)/(log_xmax-log_xmin) * ((r+1)*50) + 1 # effectively linear interpolation
end


"""
    grid_xval(n, r, idx)

Find the value of ``x`` associated with a given grid point.
"""
function grid_xval(n, r, idx)
    if idx == 0 return 0; end

	log_xmin = log10(max(10^-4, 0.01 * (sqrt(n) - sqrt(r))^2))
    log_xmax = log10(10 * (sqrt(n+1) + sqrt(r+1))^2)

	return 10^(log_xmin + (log_xmax-log_xmin)*((idx-1)/((r+1)*50)))
end


"""
    nearby_grid_xvals(n, r, idx)

Return a set of x-values surrounding the given one.
"""
function nearby_grid_xvals(n, r, idx)
	if 0 < idx < (r+1)*50+1
		return grid_xval.(n, r, (idx-1, idx, idx+1))
    elseif idx == 0
        return grid_xval.(n, r, (idx, idx+1))
    else
		return grid_xval.(n, r, (idx-1, idx))
	end
end


"""
    interpolate(x::Tuple, y::Tuple, z)

Estimate the value of ``I_{n,r}(z)`` given the provided values that surround ``z``.
"""
function interpolate(x::Tuple, y::Tuple, z)
	if length(x) == 3 # quadratic lagrange interpolation
		a = (x[1]*(y[3]-y[2]) + x[2]*(y[1]-y[3]) + x[3]*(y[2]-y[1])) / ((x[1]-x[2])*(x[1]-x[3])*(x[2]-x[3]))
		b = (y[2]-y[1])/(x[2]-x[1]) - a*(x[1]+x[2])
		c = y[1] - a*x[1]^2 - b*x[1]
		a*z^2 + b*z + c
	elseif length(x) == 2 # lerp
		(y[1]*(x[2]-z) + y[2]*(z-x[1])) / (x[2]-x[1])
	end
end


"""
    I(n, r, x)

Calculate the value of ``I_{n,r}(x)``.
"""
function I(n, r, x)
    if x <= 0 || n < 0 || r < 0
        return 0 # I(n, r, 0) = 0, and if x < 0 then x really ought to be zero
    end 
    if n < r
        return (-1)^(n-r)*I(r, n, x) # handles a case that our functions don't like, using an identity
    end

	n = Integer(n)
	r = Integer(r)
	
	x_idx = Integer(round(nearest_grid_idx(n, r, x)))
	row = find_row(n, r)
	if row > 31600
		@error "I: tried to retrieve a value for n > 250 or x > 250."
	end
	
	if x_idx == -1 return 0.0; end # x large
	
    # get grid x-values around the point in question
	itp_xvals = nearby_grid_xvals(n, r, x_idx)

	if x_idx > 1 && length(itp_xvals) == 3 # 
		itp_Ivals = (inr_table[row][x_idx-1],
					inr_table[row][x_idx], 
					inr_table[row][x_idx+1])
    elseif x_idx == 0 # lerp between (0, 0) and first point
        itp_Ivals = (0,
                    inr_table[row][x_idx+1])
	elseif x_idx == 1 # lagrange with (0, 0) and first two points
		itp_Ivals = (0,
					inr_table[row][x_idx],
					inr_table[row][x_idx+1])
	else 
		itp_Ivals = (inr_table[row][x_idx-1], 
					inr_table[row][x_idx])
	end

	interpolate(itp_xvals, itp_Ivals, x)
end

end #module I_nr