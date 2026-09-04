test_that("ContrastMultPred() handles Reverse Helmert contrasts correctly", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2"))
  )

  dat <- dat[rep(seq_len(nrow(dat)), each = 20), ]

  dat$y <- rnorm(
    nrow(dat),
    mean = with(
      dat,
      ifelse(
        A == "A1", 10,
        ifelse(A == "A2", 12, 15)
      ) +
        ifelse(B == "B2", 2, 0)
    ),
    sd = 2
  )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- stats::contr.treatment(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  out <- ContrastMultPred(
    model = mod,
    data = dat
  )

  expect_true(is.data.frame(out))

  expect_equal(
    sort(unique(out$contrast_id)),
    c(1, 2)
  )

  expect_true(
    all(c(
      "A1",
      "Mean(A2,A3)",
      "A2",
      "A3"
    ) %in% out$A)
  )

  expect_equal(
    unique(out$split_var),
    "B"
  )

  expect_equal(
    sort(unique(out$B)),
    c("B1", "B2")
  )
})

test_that("ContrastMultPred() recovers analysis data from model", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(20)
  )

  dat$y <- with(
    dat,
    10 +
      2 * (A == "A2") +
      4 * (A == "A3") +
      2 * (B == "B2") +
      2 * (A == "A2" & B == "B2")
  ) +
    stats::rnorm(
      nrow(dat),
      0,
      1
    )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  out_explicit <- ContrastMultPred(
    model = mod,
    data = dat,
    configuration = "all",
    methodological_note = FALSE
  )

  out_recovered <- ContrastMultPred(
    model = mod,
    configuration = "all",
    methodological_note = FALSE
  )

  expect_equal(
    out_recovered,
    out_explicit
  )
})
