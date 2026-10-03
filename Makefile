.PHONY: download verify-data reproduce check

download:
	bash code/download.sh

verify-data:
	sha256sum -c code/results/SHA256SUMS

reproduce: verify-data
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- reproduce

check:
	cargo fmt --manifest-path code/Cargo.toml -- --check
	cargo clippy --locked --offline --manifest-path code/Cargo.toml --all-targets -- -D warnings
	cargo test --locked --offline --manifest-path code/Cargo.toml
