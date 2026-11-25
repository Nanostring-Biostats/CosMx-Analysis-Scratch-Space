\name{spatial_model}
\alias{spatial_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Description of `spatial_model` arguments for including spatially correlated random effects in DE model. %% ~~function to do ... ~~
}
\description{
The `spatial_model` argument is an optional list argument to the `smi_de` function to allow for inclusion of a spatially correlated random effect in the DE model.
There are 2 options for 'name',
\itemize{
  \item "GP_Matern" (gaussian process, Matern covariance using spaMM package)
  \item "GP_INLA" (gaussian process with Matern prior covariance fit on a mesh using INLA package)
}
}

\usage{
# Option 1: GP_Matern
spatial_model = list(
  name = "GP_Matern",
  k_prop_n = 0.2,
  x_coord_col = "sdimx",
  y_coord_col = "sdimy",
  split_neighbors_by_colname = "Run_Tissue_name",
  spatial_random_effect = ~Matern(1 | sdimx_cluster + sdimy_cluster \%in\% Run_Tissue_name)
)

# Option 2: GP_INLA
spatial_model = list(
  name = "GP_INLA",
  x_coord_col = NULL,
  y_coord_col = NULL,
  max.edge.mesh = NULL,
  offset.mesh = c(-0.05, -0.1),
  alpha.inlaprior = 2,
  prior.range.inlaprior = NULL,
  prior.sigma.inlaprior = c(0.07, 0.5),
  constr.inlaprior = TRUE,
  spatial_random_effect = ~f(s, model = spde),
  quantiles = c(0.025, 0.5, 0.975)
)
}

\arguments{
  \item{spatial random effect follows a Gaussian Process, Matern covariance using spaMM package}{arguments when \code{name = "GP_Matern"}
    \describe{
      \item{name}{\code{"GP_Matern"}; spatial random effect follows Gaussian Process with Matern covariance structrue}
      \item{k_prop_n}{default: \code{0.2}; i.e. 0.2 * total # of cells will be the number of spatial clusters with a unique random effect}
      \item{x_coord_col}{default: \code{"sdimx"}; name of column in cell-level metadata representing x-coordinate}
      \item{y_coord_col}{default: \code{"sdimy"}; name of column in cell-level metadata representing y-coordinate}
      \item{split_neighbors_by_colname}{default: \code{"Run_Tissue_name"}; metadata column indicating flow-cell ID or sample ID}
      \item{spatial_random_effect}{default: \code{~Matern(1 | sdimx_cluster + sdimy_cluster \%in\% Run_Tissue_name)}; spatial random effect based on the x,y centroid coordinates of k-means clusters, nested within 'Run_Tissue_name' (sample ID / flow cell ID)}
    }
  }

  \item{gaussian process with Matern prior covariance fit on a mesh using INLA package}{arguments when \code{name = "GP_INLA"}
    \describe{
      \item{name}{\code{"GP_INLA"}; gaussian process with Matern prior covariance fit on a mesh using INLA package}
      \item{x_coord_col}{default: \code{NULL}; name of column in cell-level metadata representing x-coordinate. Must specify unless passing in a precomputed \code{mesh}}
      \item{y_coord_col}{default: \code{NULL}; name of column in cell-level metadata representing y-coordinate. Must specify unless passing in a precomputed \code{mesh}}
      \item{max.edge.mesh}{\code{max.edge} argument passed to \code{INLA::inla.mesh.2d}. If not specified, this is computed by default as \code{mean(c(rangex, rangey)) / sqrt(ncells) * c(5, 10)}, where \code{ncells} is the number of modeled cells, and \code{rangex}, \code{rangey} are the difference between minimum and maximum \code{x}/\code{y} coordinates, respectively.}
      \item{offset.mesh}{default: \code{c(-0.05, -0.1)}; \code{offset} argument passed to \code{INLA::inla.mesh.2d}}
      \item{alpha.inlaprior}{default: 2; \code{alpha} argument passed to \code{INLA::inla.spde2.pcmatern}}
      \item{prior.range.inlaprior}{\code{prior.range} argument passed to \code{INLA::inla.spde2.pcmatern}, If not specified, this will be computed by default as \code{c(mean(c(rangex, rangey)), 0.5)}.}
      \item{prior.sigma.inlaprior}{default: \code{c(0.07, 0.5)}; \code{prior.sigma} argument passed to \code{INLA::inla.spde2.pcmatern} }
      \item{constr.inlaprior}{default: TRUE; \code{contr} argument passed to \code{INLA::inla.spde2.pcmatern}}
      \item{spatial_random_effect}{default: \code{~f(s, model = spde)}; formula for spatial random effect used for \code{INLA::inla}}
      \item{quantiles}{default: \code{c(0.025, 0.5, 0.975)}; posterior quantiles (i.e. for credible intervals) to be computed for contrasts}
    }
  }
}

\examples{
# Option 1: GP_Matern
smi_de(
  assay_matrix
  , metadata
  , formula
  , family
  , targets = rownames(sem)[1:3]
  , spatial_model = list(
    name = "GP_Matern" # spatial random effect follows Gaussian Process with Matern covariance structure
    , k_prop_n = 0.01 # 0.01 * total # of cells will be the number of spatial clusters with a unique random effect
    , x_coord_col = "sdimx" # name of column in cell-level metadata representing x-coordinate
    , y_coord_col = "sdimy" # name of column in cell-level metadata representing y-coordinate
    , split_neighbors_by_colname = "Run_Tissue_name" # metadata column indicating flow-cell ID or sample ID
    , spatial_random_effect = ~Matern(1 | sdimx_cluster + sdimy_cluster %in% Run_Tissue_name)
  )
)

# Option 2: GP_INLA
smi_de(
  assay_matrix
  , metadata
  , formula
  , spatial_model = list(
    # default arguments for `name = "GP_INLA"`
    name = "GP_INLA"
    , x_coord_col = "sdimx"
    , y_coord_col = "sdimy"
    , max.edge.mesh = c(20000, 2000 * 3)
    , offset.mesh = c(0.2, 0.2)
    , alpha.inlaprior = 2
    , prior.range.inlaprior = c(1, 0.5)
    , prior.sigma.inlaprior = c(0.07, 0.5)
    , constr.inlaprior = TRUE
    , spatial_random_effect = ~f(s, model = spde)
    , quantiles = c(0.025, 0.5, 0.975)
  )
)
}

\details{
%% ~~ If necessary, more details than the description above ~~
}
\value{
%% ~Describe the value returned
%% If it is a LIST, use
%% \item{comp1 }{Description of 'comp1'}
%% \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%% ~~who you are~~
}
\note{
%% ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
See also:
\code{\link[spaMM]{fitme}}
\code{\link[INLA]{inla}}
\code{\link[INLA]{inla.mesh.2d}}
\code{\link[INLA]{inla.spde2.pcmatern}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
