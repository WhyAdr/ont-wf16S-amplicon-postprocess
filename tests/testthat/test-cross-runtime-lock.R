# =============================================================================
# R/Python taxonomy-cache lock interoperability
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "dependencies.R"))
source(file.path("..", "..", "analysis", "utils", "atomic_io.R"))

python <- tryCatch(find_python(), error = function(e) NA_character_)
resolver_path <- normalizePath(
  file.path("..", "..", "analysis", "utils", "ncbi_taxonomy.py"),
  winslash = "/", mustWork = TRUE
)

python_lock_script <- paste(
  "import importlib.util,sys,time",
  "spec=importlib.util.spec_from_file_location('resolver',sys.argv[1])",
  "module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)",
  "with module.acquire_cache_lock(sys.argv[2],timeout=float(sys.argv[3])):",
  " print('LOCKED',flush=True);time.sleep(float(sys.argv[4]))",
  sep = "\n"
)

test_that("an R-held taxonomy lock rejects direct Python acquisition", {
  skip_if(is.na(python), "Python is unavailable")
  root <- tempfile("r_python_lock_")
  dir.create(root)
  cache <- file.path(root, "cache.json")
  writeLines("{}", cache)
  held <- acquire_taxonomy_lock(cache, timeout_ms = 1000)
  on.exit(release_taxonomy_lock(held), add = TRUE)
  child <- processx::run(
    python, c("-c", python_lock_script, resolver_path, cache, "0.2", "0"),
    error_on_status = FALSE, timeout = 5000
  )
  expect_gt(child$status, 0L)
  expect_match(paste(child$stderr, child$stdout), "E_TAXONOMY_CACHE_BUSY")
})

test_that("a Python-held taxonomy lock rejects R acquisition", {
  skip_if(is.na(python), "Python is unavailable")
  root <- tempfile("python_r_lock_")
  dir.create(root)
  cache <- file.path(root, "cache.json")
  writeLines("{}", cache)
  child <- processx::process$new(
    python, c("-c", python_lock_script, resolver_path, cache, "1", "3"),
    stdout = "|", stderr = "|", cleanup = TRUE
  )
  on.exit(if (child$is_alive()) child$kill(), add = TRUE)
  ready <- FALSE
  for (attempt in seq_len(20L)) {
    child$poll_io(100)
    text <- paste(child$read_output_lines(), collapse = "\n")
    if (grepl("LOCKED", text, fixed = TRUE)) {
      ready <- TRUE
      break
    }
  }
  expect_true(ready)
  expect_error(acquire_taxonomy_lock(cache, timeout_ms = 100),
               "E_TAXONOMY_CACHE_BUSY")
})

test_that("cache-lock-held bypass requires an audited owner PID", {
  skip_if(is.na(python), "Python is unavailable")
  root <- tempfile("lock_bypass_")
  dir.create(root)
  abundance <- file.path(root, "abundance.tsv")
  cache <- file.path(root, "cache.json")
  writeLines(c(
    "tax\tS1",
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t1"
  ), abundance)
  writeLines("{}", cache)
  result <- processx::run(
    python, c(
      resolver_path, "--abundance", abundance, "--cache", cache,
      "--mode", "refresh", "--defer-cache-commit", "--cache-lock-held",
      "--cache-lock-owner-pid", "99999999"
    ), error_on_status = FALSE
  )
  expect_gt(result$status, 0L)
  expect_match(result$stderr, "cache-lock-owner-pid")
})
