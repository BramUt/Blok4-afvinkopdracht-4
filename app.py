from flask import Flask, render_template, request
from Bio.Seq import Seq
from pymodules import is_dna, is_prot, is_rna, get_gene

app = Flask(__name__)


@app.route('/', methods=["GET"])
def hello_world():
    """Pagina die het type van een opgegeven sequentie kan herkennen.
    Als het DNA is wordt deze naar RNA en eiwit omgezet. Is het RNA dan
    wordt het naar eiwit omgezet. Van eiwitsequenties wordt met behulp
    van blastp en biojbehoren gen gezocht.
    """
    seq = request.args.get("seq", "").upper()
    dna_seq = ""
    rna_seq = ""
    prot_seq = ""
    # Kijkt of de sequentie DNA is.
    if is_dna(seq):
        seq_type = "DNA"
        dna_object = Seq(seq)
        # Transcriptie van sequentie
        rna_seq = Seq.transcribe(dna_object)
        # Translatie van sequentie.
        prot_seq = Seq.translate(rna_seq)
        return render_template("afvink4page.html", seq=seq, seq_type=seq_type,
                               dna_seq=seq, rna_seq=str(rna_seq),
                               prot_seq=str(prot_seq))
    # Kijkt of de sequentie RNA is.
    elif is_rna(seq):
        seq_type = "RNA"
        seq_object = Seq(seq)
        # Translatie van de sequentie.
        prot_seq = Seq.translate(seq_object)
        return render_template("afvink4page.html", seq=seq, seq_type=seq_type,
                               rna_seq=seq, prot_seq=str(prot_seq))
    # Kijkt of de sequenite een eiqit is.
    elif is_prot(seq):
        seq_type = "Protein"
        # Gen wordt opgezocht.
        # TODO: Niet kunnen test, BLAST blijft oneindig door gaan.
        accessie_code = get_gene(seq)
        return render_template("afvink4page.html", seq=seq, seq_type=seq_type,
                               prot_seq=seq, accessie_code=accessie_code)
    # Wordt aangeroepen als de sequentie niet herkend wordt.
    else:
        seq_type = "Unrecognized"
        return render_template("afvink4page.html", seq=seq, seq_type=seq_type)


if __name__ == '__main__':
    app.run()
