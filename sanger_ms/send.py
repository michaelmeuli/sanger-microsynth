"""Send the PDF report by SMTP (STARTTLS)."""
from __future__ import annotations

import smtplib
from email.message import EmailMessage
from pathlib import Path


def send_report(cfg: dict[str, str], pdf: Path, subject: str, body: str, extra: list[Path] = ()) -> None:
    msg = EmailMessage()
    msg["From"], msg["To"], msg["Subject"] = cfg["MAIL_FROM"], cfg["MAIL_TO"], subject
    msg.set_content(body)
    for f in [pdf, *extra]:
        msg.add_attachment(f.read_bytes(), maintype="application", subtype="pdf", filename=f.name)
    with smtplib.SMTP(cfg["SMTP_HOST"], int(cfg["SMTP_PORT"])) as smtp:
        smtp.starttls()
        smtp.login(cfg["SMTP_USER"], cfg["SMTP_PASSWORD"])
        smtp.send_message(msg)
