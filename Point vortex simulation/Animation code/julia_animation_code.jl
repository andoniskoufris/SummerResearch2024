using Plots
using FileIO
using DelimitedFiles
using LaTeXStrings


function read_particle_data(filename)
    data = readdlm(filename)
    num_particles = size(data, 2) ÷ 2
    return reshape(data, :, num_particles, 2)
end


function create_animation(particle_data, output_folder, output_file)
    num_frames = floor(Int, size(particle_data, 2)/2)
    num_particles = size(particle_data, 1)
    radius = 1.05 * maximum([ maximum(sqrt.(particle_data[:,1].^2 + particle_data[:,2].^2)), maximum(sqrt.(particle_data[:,num_frames-1].^2 + particle_data[:,num_frames].^2)) ] )
    
    animation = @animate for i in 1:2:floor(Int, num_frames)
        plot(particle_data[:, i], particle_data[:, i+1], seriestype=:scatter, legend=false, xlimits=(-radius, radius), ylimits=(-radius, radius), aspect_ratio=:equal, color="royalblue1", markersize=1, dpi=300)
        #xlabel!(L"$x-\mathrm{axis}$")
        #ylabel!(L"$y-\mathrm{axis}$")
    end

    # Combine the output folder and file name to get the full path
    full_output_path = joinpath(output_folder, output_file)
    
    # Save the animation as an MP4 file in the specified folder
    mp4(animation, full_output_path, fps=30)

end

NumVortices = 4027

# Specify the full path to your data file
data_file = raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\\" * string(NumVortices) * "_single_perturbation_1.2.txt"

# Read particle data from the file
particle_data = readdlm(data_file)

# Specify the output file name
output_file = string(NumVortices) * raw"_test.mp4"

# Specify the full path to the output folder
output_folder = raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Animations"

# Create and save the animation as an MP4 file in the specified folder
create_animation(particle_data, output_folder, output_file)