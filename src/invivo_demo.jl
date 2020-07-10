const gbm_study_uid = "1.3.6.1.4.1.14519.5.2.1.4591.4001.304604545029494418165835320551"

function download_invivo_studies(study=gbm_study_uid; destination::AbstractString, overwrite = false)
    make_folder(destination)
    (vfa_folder, dce_folder) = download_vfa_and_dce(study; destination, overwrite)
    return (vfa_folder, dce_folder)
end

function download_vfa_and_dce(study_id; destination, overwrite = false)
    vfa_folder = joinpath(destination, study_id, "vfa")
    dce_folder = joinpath(destination, study_id, "dce")

    if isdir(vfa_folder) && isdir(dce_folder)
        if !overwrite
            return (vfa_folder, dce_folder)
        end
    end

    gbm_series = tcia_series(study = study_id)
    vfa_series = find_vfa_series(gbm_series)
    dce_series = find_dce_series(gbm_series)

    if !isdir(vfa_folder) || overwrite == true
        download_series(vfa_series, destination = vfa_folder)
    end

    if !isdir(dce_folder) || overwrite == true
        download_series(dce_series, destination = dce_folder)
    end
    return (vfa_folder, dce_folder)
end

find_vfa_series(series_dataframe) = find_in_description("MAP", series_dataframe)
find_dce_series(series_dataframe) = find_in_description("DYN", series_dataframe)

function find_in_description(word_to_find::AbstractString, series_dataframe)
    descriptions = series_dataframe.SeriesDescription
    found_indices = findall(occursin.(word_to_find, descriptions))
    if length(found_indices) < 1
        errant_id = series_dataframe.PatientID[1]
        errant_study = series_dataframe.StudyInstanceUID[1]
        error("No single series in $errant_id found with $word_to_find in their description.
              The full study UID is $errant_study\n")
    elseif length(found_indices) > 1
        errant_id = series_dataframe.PatientID[1]
        @warn "Multiple series in $errant_id found containing $word_to_find in their description.
              An arbitrary series will be used among the found series"
        found_index = found_indices[end]
    else
        found_index = found_indices[1]
    end
    found_series = series_dataframe.SeriesInstanceUID[found_index]
    return found_series
end

function download_series(series_id::AbstractString; destination)
    make_folder(destination; remove_existing = true)
    zip_file = joinpath(destination, "downloaded.zip")
    tcia_images(series = series_id, file = zip_file)
    unzip_command = `unzip -o $zip_file -d $destination`
    run(unzip_command)
    rm(zip_file)
    return destination
end

function make_folder(desired_folder::AbstractString; remove_existing = false)
    if !isdir(desired_folder)
        mkpath(desired_folder)
    elseif remove_existing
        rm(desired_folder, recursive = true)
        mkpath(desired_folder)
    end
    return desired_folder
end