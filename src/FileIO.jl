"""
    writeMps!

Write the deterministic equivalent formulation in MPS format.
"""
function writeMps!(m::SJ.StructuredModel, filename::AbstractString; options...)
    # free any existing model pointer
    freeModel(dspenv)

    # set options
    setoptions!(options)

    # classify the second stage constraints into linear ones and quadratic ones
    classifyConstrs(m)

    # load problem
    load_problem!(m)

    # write problem into MPS file
    writeMps!(dspenv, filename)

    return
end