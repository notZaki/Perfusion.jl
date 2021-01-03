function load_vfa_dicom(dir)
    dicoms = load_dicom_dir(dir)
    num_images = length(dicoms)
    unique_flip_angles = unique(lookup.(dicoms, "Flip Angle"))
    num_flip_angles = length(unique_flip_angles)
    num_slices = Int(num_images / num_flip_angles)

    dummy_image = lookup(dicoms[1], "Pixel Data")
    image_size = size(dummy_image)
    signal = zeros(num_images, image_size...)
    flip_angles = zeros(num_images)
    for dicom in dicoms
        instance = lookup(dicom, "Instance Number")
        signal[:, :, instance] = lookup(dicom, "Pixel Data")
        flip_angles[instance] = lookup(dicom, "Flip Angle")
    end
    signal = reshape(signal, (image_size..., num_slices, num_flip_angles))
    flip_angles = reshape(flip_angles, (num_slices, num_flip_angles))[1, :]
    @. flip_angles = deg2rad(flip_angles)
    TR = lookup(dicoms[1], "Repetition Time")
    return (; signal, angles = flip_angles, TR)
end

function load_dce_dicom(dir; num_slices::Integer)
    dicom_files = load_dicom_dir(dir)
    num_images = length(dicom_files)
    num_timepoints = Int(num_images / num_slices)

    dummy_image = lookup(dicom_files[1], "Pixel Data")
    image_size = size(dummy_image)
    signal = zeros(image_size..., num_images)
    timepoints = zeros(num_images)
    for dicom in dicom_files
        instance = lookup(dicom, "Instance Number")
        signal[:, :, instance] = lookup(dicom, "Pixel Data")
        timepoints[instance] = lookup(dicom, "Trigger Time")
    end
    signal = reshape(signal, (image_size..., num_slices, num_timepoints))
    timepoints = reshape(timepoints, (num_slices, num_timepoints))[1, :] ./ 1000 ./ 60
    TR = lookup(dicom_files[1], "Repetition Time")
    flip_angle = deg2rad(lookup(dicom_files[1], "Flip Angle"))
    return (; signal, t = vec(timepoints), TR, angle = flip_angle)
end

function load_dicom_dir(dir)
    dicom_files = list_dicom_in_dir(dir)
    num_files = length(dicom_files)
    dummy_data = dcm_parse(dicom_files[1])
    dicom_data = Vector{typeof(dummy_data)}(undef, num_files)
    for (index, file) in enumerate(dicom_files)
        dicom_data[index] = dcm_parse(file)
    end
    return dicom_data
end

function list_dicom_in_dir(dir)
    files = readdir(dir)
    dicom_files = joinpath.(dir, files[is_dicom.(files)])
    return dicom_files
end

is_dicom(file::AbstractString) = splitext(file)[2] == ".dcm"
