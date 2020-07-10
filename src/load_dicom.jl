function load_vfa_dicom(; folder)
    dicom_data = load_dicom_from_folder(folder)
    number_of_images = length(dicom_data)
    unique_flip_angles = unique(lookup.(dicom_data, "Flip Angle"))
    number_of_flip_angles = length(unique_flip_angles)
    number_of_slices = Int(number_of_images / number_of_flip_angles)

    dummy_image = lookup(dicom_data[1], "Pixel Data")
    image_size = size(dummy_image)
    signal_data = zeros(image_size..., number_of_images)
    flip_angles = zeros(number_of_images)
    for dicom in dicom_data
        instance = lookup(dicom, "Instance Number")
        signal_data[:,:,instance] = lookup(dicom, "Pixel Data")
        flip_angles[instance] = lookup(dicom, "Flip Angle")
    end
    signal_data = reshape(signal_data, (image_size..., number_of_slices, number_of_flip_angles))
    flip_angles = reshape(flip_angles, (number_of_slices, number_of_flip_angles))[1,:]
    @. flip_angles = deg2rad(flip_angles)
    TR = lookup(dicom_data[1], "Repetition Time")
    return (signal = signal_data, angles = flip_angles, TR = TR)
end

function load_dce_dicom(; folder, num_slices::Int = 16)
    dicom_data = load_dicom_from_folder(folder)
    number_of_images = length(dicom_data)
    number_of_timepoints = Int(number_of_images / num_slices)

    dummy_image = lookup(dicom_data[1], "Pixel Data")
    image_size = size(dummy_image)
    signal_data = zeros(image_size..., number_of_images)
    timepoints = zeros(number_of_images)
    for dicom in dicom_data
        instance = lookup(dicom, "Instance Number")
        signal_data[:,:,instance] = lookup(dicom, "Pixel Data")
        timepoints[instance] = lookup(dicom, "Trigger Time")
    end
    signal_data = reshape(signal_data, (image_size..., num_slices, number_of_timepoints))
    timepoints = reshape(timepoints, (num_slices, number_of_timepoints))[1,:] ./ 1000 ./ 60
    TR = lookup(dicom_data[1], "Repetition Time")
    flip_angle = deg2rad(lookup(dicom_data[1], "Flip Angle"))
    return (signal = signal_data, timepoints = vec(timepoints), TR = TR, angle = flip_angle)
end

function load_dicom_from_folder(folder)
    dicom_files = get_dicom_files_in_folder(folder = folder)
    number_of_files = length(dicom_files)
    dummy_data = dcm_parse(dicom_files[1])
    dicom_data = Vector{typeof(dummy_data)}(undef, number_of_files)
    for (index, file) in enumerate(dicom_files)
        dicom_data[index] = dcm_parse(file)
    end
    return dicom_data
end

function get_dicom_files_in_folder(; folder)
    files = readdir(folder)
    dicom_files = joinpath.(folder, files[is_dicom.(files)])
    return dicom_files
end

is_dicom(file::AbstractString) = splitext(file)[2] == ".dcm"