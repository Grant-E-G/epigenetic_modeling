.PHONY: download verify-data reproduce followup check

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

# Extended validation uses generated paired subsets and full raw matrices.
followup: verify-data
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- arrays
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- forecast
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- kinetic-followup
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- growth-followup
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- matching-followup
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- followup-report
