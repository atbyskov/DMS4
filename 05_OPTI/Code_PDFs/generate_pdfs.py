import os
from fpdf import FPDF

BASE = r"C:\Users\boggy\OneDrive\Desktop\DMS4\05_OPTI"
OUT = os.path.join(BASE, "Code_PDFs")

FILES = [
    "Main.py",
    "MyAPDLCall.py",
    "Post_Process.py",
    "optimization.py",
    "APDL_Input.py",
    "aggregate.py",
    "acs.py",
]

class CodePDF(FPDF):
    def __init__(self, title):
        super().__init__(orientation="P", unit="mm", format="A4")
        self._title = title

    def header(self):
        self.set_font("Courier", "B", 10)
        self.cell(0, 8, self._title, border=0, align="C", new_x="LMARGIN", new_y="NEXT")
        self.ln(2)

    def footer(self):
        self.set_y(-15)
        self.set_font("Courier", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="C")


for fname in FILES:
    src = os.path.join(BASE, fname)
    with open(src, "r", encoding="utf-8") as f:
        lines = f.readlines()

    pdf = CodePDF(title=fname)
    pdf.alias_nb_pages()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)

    num_width = 12
    line_height = 4.2
    code_width = 190 - num_width - 2

    for i, line in enumerate(lines, 1):
        text = line.rstrip("\n\r")
        text = text.replace("\t", "    ")
        text = text.encode("latin-1", errors="replace").decode("latin-1")

        pdf.set_font("Courier", "B", 8)
        pdf.set_text_color(120, 120, 120)
        pdf.cell(num_width, line_height, f"{i:>4}", border=0)

        pdf.set_font("Courier", "", 8)
        pdf.set_text_color(0, 0, 0)
        pdf.cell(2, line_height, " ", border=0)

        if len(text) > 115:
            chunks = [text[j:j+115] for j in range(0, len(text), 115)]
            for ci, chunk in enumerate(chunks):
                if ci > 0:
                    pdf.cell(num_width + 2, line_height, "", border=0)
                pdf.cell(code_width, line_height, chunk, new_x="LMARGIN", new_y="NEXT")
        else:
            pdf.cell(code_width, line_height, text, new_x="LMARGIN", new_y="NEXT")

    out_name = os.path.splitext(fname)[0] + ".pdf"
    out_path = os.path.join(OUT, out_name)
    pdf.output(out_path)
    print(f"Created: {out_path}")

print("\nDone! All PDFs saved to:", OUT)
