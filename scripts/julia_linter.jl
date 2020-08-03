#!/usr/bin/env julia
using Lint
using LightXML

#=
    Code check prints xml checkStyle
=#

function main(args)

    xdoc = XMLDocument()
    xroot = create_root(xdoc, "checkstyle")

    for file in args
        file_root = new_child(xroot, "file")

        errs = collect(lintfile(file))
        if length(errs) > 0
            set_attribute(file_root, "name", file)
            for err in errs
                xs1 = new_child(file_root, "error")
                xml_info = Dict("line"=>err.line, "column" => "1",
                            "severity"=>"warning",
                            "message"=>err.scope*"(): "*err.variable*": "*err.message*"("*string(err.code)*")",
                            "source"=>"com.puppycrawl.tools.checkstyle.checks.AvoidEscapedUnicodeCharactersCheck")
                set_attributes(xs1, xml_info)
            end
        end
    end
    save_file(xdoc, "julia_lint-sonar.xml")

end

main(ARGS)
