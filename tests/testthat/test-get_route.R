# get route
data(liver_1.3_network)
data(liver_1.3_rwr_closest_dfr)
data(signature_maison)

signature_vids <- signature_maison$acetaminophen_all_all
res_diffusion <- get_route(liver_1.3_network, liver_1.3_rwr_closest_dfr, signature_vids)
res_report <- report(res_diffusion)

test_that("get_route works", {
    expect_s3_class(res_diffusion, "get_route.res")
})

test_that("report works", {
    expect_type(res_report, "list")
    expect_equal(length(res_report), 9)
})
