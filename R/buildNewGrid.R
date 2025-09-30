#' Build a new grid given two coordinate matrices
#'
#' Given two coordinate matrices, concave hulls are found around them. The intersection
#' between these two hulls is found, and in that area an evenly spaced, discrete
#' grid is constructed.
#'
#' This function is mainly used to create a grid on which fitted GAMs can be
#' evaluated to calculate correlations.
#'
#' @param Cx,Ey The coordinate matrices
#' @param n_points_grid an integer, the number of points desired for the new grid
#' @details The new grid will contain approximately the number of new points requested,
#' depending on the size of the concave hull
#' @returns A data frame of two columns with all points of the grid,
#' with column names x and y.
#' @export
#' @importFrom concaveman concaveman
#' @importFrom sf st_area st_sf st_sfc st_polygon st_intersection st_make_grid st_within st_coordinates
#' @seealso \link[concaveman]{concaveman}
#' @examples
#' Cx = matrix(runif(40, 0, 1), 20, 2)
#' Ey = matrix(runif(30, 0, 1), 15, 2)
#' buildNewGrid(Cx, Ey, n_points_grid = 50)
buildNewGrid = function(Cx, Ey, n_points_grid){
    stopifnot(ncol(Cx) == 2L, ncol(Ey)==2, is.numeric(n_points_grid))
    #Find concave hulls
    hullx <- concaveman(Cx)
    hully <- concaveman(Ey)
    # Create sf polygons
    polygon_sf_x <- st_sf(st_sfc(st_polygon(list(hullx))))
    polygon_sf_y <- st_sf(st_sfc(st_polygon(list(hully))))
    polygon_sf = st_intersection(polygon_sf_x, polygon_sf_y)
    # Estimate grid spacing based on area and target point count
    area <- st_area(polygon_sf)
    grid_spacing <- sqrt(as.numeric(area) / n_points_grid)
    # Create a regular grid within the bounding box
    grid <- st_make_grid(polygon_sf, what = "centers", cellsize = grid_spacing)
    # Keep only points that fall inside the polygon
    grid_inside <- grid[st_within(grid, polygon_sf, sparse = FALSE)[,1]]
    grid_inside <- st_coordinates(grid_inside)[, c("X", "Y")]
    colnames(grid_inside) = c("x", "y")
    return(as.data.frame(grid_inside))
}
