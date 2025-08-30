using Plots
using FileIO
using DelimitedFiles
using LaTeXStrings

function read_particle_data(filename)
    data = readdlm(filename)
    num_particles = size(data, 2) ÷ 2
    return reshape(data, :, num_particles, 2)
end

function create_animation(real_coefficients, fourier_coefficients, vortex_positions, output_folder, output_file, NumEdgeVortices)
    num_frames = floor(Int, size(real_coefficients, 1)/2)
    num_particles = size(real_coefficients, 2)
    radius = 1.1 * maximum([ maximum(sqrt.(vortex_positions[1:num_particles-1,1].^2 + vortex_positions[1:num_particles-1,2].^2)), maximum(sqrt.(vortex_positions[1:num_particles-1,num_frames-1].^2 + vortex_positions[1:num_particles-1,num_frames].^2)) ] )

    max_radius_plot = maximum(real_coefficients[2:2:2*num_frames,:])
    min_radius_plot = minimum(real_coefficients[2:2:2*num_frames,:])

    max_fourier_plot = maximum(fourier_coefficients[2:2:2*num_frames,:])
    min_fourier_plot = minimum(fourier_coefficients[2:2:2*num_frames,:])

    animation = @animate for i in 1:2:num_frames

        l = @layout [[grid(2,1)] a{0.9h}]

        p1 = plot(real_coefficients[i,:],  real_coefficients[2*i,:], xlabel=L"$\theta$", ylabel=L"$R-R_0\ \mathrm{(arb. units)}$", legend=false, ylimits=(min_radius_plot, max_radius_plot))
        p2 = plot(fourier_coefficients[i,:],  fourier_coefficients[2*i,:], xlabel=L"$k\ \mathrm{(arb. units)}$", ylabel=L"\mathrm{Fourier\ coefficient}", legend=false, ylimits=(min_fourier_plot, max_fourier_plot))

        #p3 = plot(vortex_positions[:, i], vortex_positions[:, i+1], seriestype=:scatter, legend=false, xlimits=(-radius, radius), ylimits=(-radius, radius), aspect_ratio=:equal, color="royalblue1", markersize=1.5)
        p3 = plot(vortex_positions[1:NumEdgeVortices, i], vortex_positions[1:NumEdgeVortices, i+1], seriestype=:scatter, legend=false, xlimits=(-radius, radius), ylimits=(-radius, radius), aspect_ratio=:equal, color="firebrick3", markersize=1)

        plot(p1, p2, p3, layout = l, size="(1920, 1080)")
        #xlabel!(L"$x-\mathrm{axis}$")
        #ylabel!(L"$y-\mathrm{axis}$")
    end

    # Combine the output folder and file name to get the full path
    full_output_path = joinpath(output_folder, output_file)

    # Save the animation as an MP4 file in the specified folder
    mp4(animation, full_output_path, fps=15)

end

NumVortices = 4027

# Specify the full path to your data file
data_file_real = raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\\" * string(NumVortices) * "_real_coefficient.txt"
data_file_fourier = raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\\" * string(NumVortices) * "_fourier_coefficient.txt"
data_file_positions = raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\\" * string(NumVortices) * "_single_perturbation_1.2.txt"

# Read particle data from the file
real_coefficients = readdlm(data_file_real)
fourier_coefficients = readdlm(data_file_fourier)
vortex_positions = readdlm(data_file_positions)

# Specify the output file name
output_file = string(NumVortices) * raw"_edge_modes.mp4"

# Specify the full path to the output folder
output_folder = raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Animations"

# Create and save the animation as an MP4 file in the specified folder
create_animation(real_coefficients, fourier_coefficients, vortex_positions, output_folder, output_file, 218)