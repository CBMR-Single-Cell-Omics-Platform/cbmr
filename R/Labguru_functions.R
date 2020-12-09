#' Download metadata from Labguru
#'
#' Download metadata from Labguru when given a experiment id or
#' an experiment name
#'
#' @param exp_id numeric experimental id, fastest way to get data
#' @param exp_name character of experimental name
#' @export
#' @examples
#' \dontrun{
#' labguru_set_token(
#'   token = "abcdefg",
#'   server = "https://sund.labguru.com"
#' )
#' get_metadata(exp_id = 1)
#' }
#'
get_metadata <- function(exp_id = NULL, exp_name = NULL) {
  if (is.null(exp_id)) {
    if (is.null(exp_name)) stop("exp_id or exp_name must be set.")
    all_exps <- LabguruR::labguru_list_experiments()
    exp_id <- all_exps[all_exps[["title"]] == exp_name, "id"]
    if (length(exp_id) == 0) {
      stop("exp_name matched no experiments. Check name is correct, or supply exp_id instead")
    }
    if (length(exp_id) > 1) {
      stop("exp_name matched multiple experiments. Supply exp_id instead")
    }
  }

  procedures <- LabguruR::labguru_get_experiment(exp_id)$experiment_procedures
  procedure_id <- procedures[
    procedures[["experiment_procedure.name"]] == "Procedure",
    "experiment_procedure.id",
    drop = TRUE
  ]

  # this function is not exported from the labguru package, but the package has
  # not been updated to match the current API. When the package updates
  # labguru_get_experiment_procedure() will replace this function
  labguru_get_by_id <- function(type,
                                id,
                                server = Sys.getenv("LABGURU_SERVER"),
                                token = Sys.getenv("LABGURU_TOKEN")) {
    base_url <- server
    path <- paste0("/api/v1/", type, "/", id)
    query <- paste0("token=", token)
    url <- httr::modify_url(url = base_url, path = path, query = query)
    url <- utils::URLencode(url)
    resp <- httr::GET(url)
    if (httr::http_type(resp) != "application/json") {
      stop("API did not return JSON", call. = FALSE)
    }
    parsed <- jsonlite::fromJSON(httr::content(resp, as = "text"),
      simplifyVector = FALSE,
      simplifyDataFrame = TRUE,
      flatten = TRUE
    )
    if (httr::http_error(resp)) {
      stop(sprintf(
        "API request failed [%s]\n%s", parsed$status,
        parsed$error
      ), call. = FALSE)
    }
    parsed
  }

  parsed <- labguru_get_by_id(
    type = "sections",
    id = procedure_id
  )$elements
  sheet_id <- parsed[parsed[["element_type"]] == "excel", "id", drop = TRUE]

  datasheet_JSON <- LabguruR::labguru_get_element(sheet_id)$data

  datasheet_list <- jsonlite::fromJSON(
    jsonlite::fromJSON(txt = datasheet_JSON)$spread
  )
  datasheet_list <- datasheet_list$sheets$Sheet1$data$dataTable
  get_values <- function(x) {
    lapply(x, magrittr::extract2, "value")
  }
  filter_fn <- function(x) !is.null(x$`1`)

  datasheet <- lapply(datasheet_list, get_values)
  datasheet <- Filter(f = filter_fn, datasheet)

  out <- as.data.frame(matrix(
    nrow = length(datasheet) - 1,
    ncol = length(datasheet$`0`)
  ))
  colnames(out) <- unlist(datasheet$`0`)
  datasheet <- datasheet[-1]

  for (i in seq_along(datasheet)) {
    for (j in names(datasheet[[i]])) {
      idx <- as.integer(j) + 1
      val <- datasheet[[i]][[j]]
      if (is.null(val)) val <- NA
      out[i, idx] <- val
    }
  }
  data.table::setDT(out)
  
  out[]
}
