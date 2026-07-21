module GridExportVTK

using MultiComponentReactiveMixtureProject, ExtendableGrids

function ExportGrid(;dim=3,nref=0)
    grid, inlet_boundaries, irradiated_boundaries,
     outlet_boundaries, side_boundaries, catalyst_regions = PTR_grid_boundaries_regions(dim,nref=nref)    

    nn = num_nodes(grid)
    nbf = num_bfaces(grid)

    # Boundary Face Nodes, save as point data to VTK
    point_data_region_1 = zeros(nn)

    # region 1
    reg = 1
    for n in 1:nn
        for (i,bf) in enumerate(eachcol(grid[BFaceNodes]))
            if grid[BFaceRegions][i] == reg # 
                if n in bf
                    point_data_region_1[n] = reg
                end
            end
        end
    end
    #return point_data_region_1

    ExtendableGrids.writeVTK("out/PTR_grid3D_regions_boundaries", grid, 
        cell_data = grid[CellRegions],
        region_1 = point_data_region_1
    )
    return nothing
end

end