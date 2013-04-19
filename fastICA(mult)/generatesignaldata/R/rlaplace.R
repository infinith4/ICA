rlaplace <-
function (n, location = 0, scale = 1) 
{
    if (!is.Numeric(n, positive = TRUE, integer.valued = TRUE, 
        allowable.length = 1)) 
        stop("bad input for argument 'n'")
    if (!is.Numeric(scale, positive = TRUE)) 
        stop("'scale' must be positive")
    location = rep(location, length.out = n)
    scale = rep(scale, length.out = n)
    r = runif(n)
    location - sign(r - 0.5) * scale * log(2 * ifelse(r < 0.5, 
        r, 1 - r))
}
