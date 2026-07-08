# Seed copy of the Homebrew formula. The canonical copy lives in
# bpp/homebrew-tap (Formula/bpp-seqs.rb); add this there once, and the
# release workflow's bump-homebrew job keeps its url + sha256 current on
# every tagged release. The sha256 below is a placeholder for the initial
# seed -- the first `v*` tag push overwrites it with the real checksum.
class BppSeqs < Formula
  desc "Convert sequence/variant data (BAM/CRAM, gVCF, FASTA, PHYLIP, NEXUS) to BPP format"
  homepage "https://github.com/bpp/bpp-seqs"
  url "https://github.com/bpp/bpp-seqs/archive/refs/tags/v0.1.0.tar.gz"
  sha256 "0000000000000000000000000000000000000000000000000000000000000000"
  license "AGPL-3.0-or-later"
  head "https://github.com/bpp/bpp-seqs.git", branch: "main"

  depends_on "htslib"

  def install
    system "make"
    bin.install "bpp-seqs"
    pkgshare.install "examples"
    doc.install "README.md"
  end

  test do
    assert_match "bpp-seqs", shell_output("#{bin}/bpp-seqs --version")
    # The FASTA example inspects cleanly (no external tools required).
    ex = pkgshare/"examples/01-fasta"
    system bin/"bpp-seqs", "--dry-run", ex/"locusA.fa", ex/"locusB.fa",
           "--imap", ex/"imap.txt"
  end
end
