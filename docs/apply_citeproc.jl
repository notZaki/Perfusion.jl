function find_files_by_ext(startdir, ext = [""])
    foundfiles = String[]
    for (root, dirs, files) in walkdir(startdir)
        for file in files
            if last(splitext(file)) ∈ ext
                push!(foundfiles, joinpath(root, file))
            end
        end
    end
    return normpath.(foundfiles)
end

find_files_by_ext(startdir, ext::AbstractString) = find_files_by_ext(startdir, [ext])

function apply_citeproc(startdir)
    bib = first(find_files_by_ext(startdir, ".bib"))
    csl = first(find_files_by_ext(startdir, ".csl"))
    html = find_files_by_ext(startdir, ".html")
    for file in html
        run(`pandoc $file -o $file --bibliography $bib --csl $csl --metadata link-citations=true --filter pandoc-citeproc -f markdown --to html `)
    end
    return
end
