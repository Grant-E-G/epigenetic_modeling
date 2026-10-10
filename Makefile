.PHONY: download verify-data reproduce followup validation recovery broad-download broad-validation power dog-download dog-validation check

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

# Conversion tools live under ignored data/tools; install per code/README.md.
validation: verify-data
	data/tools/format-env/bin/python code/convert_public_data.py
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- validation data/raw code/results all

# Recovery pilot: uses the frozen human feature table produced by reproduce/annotate.
recovery: verify-data
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- validation data/raw code/results recovery

broad-download:
	bash code/download.sh --broad

broad-validation:
	data/tools/format-env/bin/python code/convert_public_data.py broad-verify
	data/tools/format-env/bin/python code/convert_public_data.py broad
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- validation data/raw code/results broad

# Rebuild the exact observed audit design before conditional resampling.
power: broad-validation
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- validation data/raw code/results power

# Frozen dog cohort counts/reference files; format environment is described in code/README.md.
dog-download:
	data/tools/format-env/bin/python code/dog_aging_data.py download

dog-validation:
	data/tools/format-env/bin/python code/dog_aging_data.py decode
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- dog-aging data/raw code/results all
	cargo run --release --locked --offline --manifest-path code/Cargo.toml -- dog-aging-audit data/raw code/results
