package(default_visibility = ["//test:__subpackages__"])

cc_binary(
    name = "SVDQuest",
    srcs = ["SVDQuest.cpp"],
    deps = [":SVDQuestTripartitionScorer",
        "@phylonaut",
        "@phylokit",
    	"@com_github_google_glog//:glog",
    ]
)

cc_library(
    name = "SVDQuestTripartitionScorer",
    hdrs = ["SVDQuestTripartitionScorer.hpp"],
    srcs = ["SVDQuestTripartitionScorer.cpp"],
    deps = ["@phylonaut", "@phylokit"]
)