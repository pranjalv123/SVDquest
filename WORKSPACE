load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

git_repository(
    name = "com_github_nelhage_rules_boost",
    commit = "9f9fb8b2f0213989247c9d5c0e814a8451d18d7f",
    remote = "https://github.com/nelhage/rules_boost",
    shallow_since = "1570056263 -0700",
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")

boost_deps()

http_archive(
    name = "com_github_gflags_gflags",
    sha256 = "34af2f15cf7367513b352bdcd2493ab14ce43692d2dcd9dfc499492966c64dcf",
    strip_prefix = "gflags-2.2.2",
    urls = ["https://github.com/gflags/gflags/archive/v2.2.2.tar.gz"],
)

http_archive(
    name = "com_github_google_glog",
    sha256 = "62efeb57ff70db9ea2129a16d0f908941e355d09d6d83c9f7b18557c0a7ab59e",
    strip_prefix = "glog-d516278b1cd33cd148e8989aec488b6049a4ca0b",
    urls = ["https://github.com/google/glog/archive/d516278b1cd33cd148e8989aec488b6049a4ca0b.zip"],
)

git_repository(
    name = "catch2",
    commit = "c5538476052dfe9d3ff2325198b1a8f42fc10669",
    remote = "https://github.com/evanmoran/catch2-bazel/",
    shallow_since = "1530732979 -0700",
)

git_repository(
    name = "phylokit",
    commit = "aac94d561d24cb747a1e8381ffa2da0a0a0d72b6",
    remote = "https://github.com/pranjalv123/phylokit/",
)

git_repository(
    name = "phylonaut",
    commit = "c49017a05a55380c6e0b04c789074bdd997eaccd",
    remote = "https://github.com/pranjalv123/phylonaut/",
)