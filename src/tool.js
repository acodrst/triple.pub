const emoji = "☘️";
const domain = "triple.pub";
const backup = Deno.env.get("CL_TRI_BACKUP");
const dt = new Date();
const tss = dt.toISOString().replaceAll(":", "").replaceAll("-", "").replaceAll(
  ".",
  "",
);
import * as base64 from "byte-base64";
import { fpng } from "fpng";
const site = { f: {} };
const files = {
  "abstract.html": "",
  "legal.html": "",
  "logo.html": "",
  "mininav.html": "",
  "pageops.html": "",
  "page.html": "",
  "pdf.html": "",
  "style.css": "",
  "subtitle.html": "",
  "thanks.html": "",
  "titlehead.html": "",
  "toc.html": "",
  "topgraph.html": "",
};
for (const file in files) {
  site.f[file] = Deno.readTextFileSync(`assets/${file}`);
}
for (let i = 0; i < 40; i++) {
  site.f[`${i}.txt`] = Deno.readTextFileSync(`assets/${i}.txt`);
  site.f[`${i}.html`] = Deno.readTextFileSync(`assets/${i}.html`);
}
site.slugs = {};
for (let i = 0; i < 40; i++) {
  const ns = parseInt(i);
  site.slugs["page" + ns] = {};
  site.slugs["page" + ns].content = Deno.readTextFileSync(
    `assets/${ns}.html`,
    "utf8",
  );
  site.slugs["page" + ns].title = Deno.readTextFileSync(
    `assets/${ns}.txt`,
    "utf8",
  ).trim();
  site.slugs["page" + ns].slug = site.slugs["page" + ns].title.normalize("NFD")
    .toLowerCase().replace(/[^a-z0-9 ]/g, "").replace(/\s+/g, "-");
}

function arr_to_hex(u8arr) {
  return `${
    Array.from(u8arr, (i) => i.toString(16).padStart(2, "0")).join("")
  }`;
}

Deno.writeTextFileSync("site.txt", `let site=${JSON.stringify(site)}\n`);
const text = Deno.readTextFileSync("site.txt") +
  Deno.readTextFileSync("dist/app.bundle.js");

const last_hash = Deno.readTextFileSync("data_sha512.txt");
const cur_hash = arr_to_hex(
  new Uint8Array(
    await crypto.subtle.digest("SHA-512", new TextEncoder().encode(text)),
  ),
);

if (last_hash.trim() != cur_hash.trim()) {
  Deno.writeTextFileSync("data_sha512.txt", cur_hash);
  const fp_obj=fpng(` Verify sig at floppypng.com - ${tss}`,text)
  const a32h = arr_to_hex(fp_obj.im.slice(-20, -16));
  console.log(`Generated FloppyPNG Size=${fp_obj.ln}`);

  const priv = Deno.readTextFileSync(Deno.env.get("CL_PRIV")).replace(
    /.*KEY-----(.+?)-----END.*/smg,
    "$1",
  );
  const b_der_str = globalThis.atob(priv);
  const b_der = Uint8Array.from([...b_der_str].map((c) =>
    c.charCodeAt()
  )).buffer;
  const prv = await globalThis.crypto.subtle.importKey(
    "pkcs8",
    b_der,
    {
      name: "RSA-PSS",
      hash: "SHA-256",
    },
    true,
    ["sign"],
  );
  const sig = await crypto.subtle.sign(
    {
      name: "RSA-PSS",
      hash: "SHA-256",
      saltLength: 32,
    },
    prv,
    fp_obj.im,
  );
  const u8sig = new Uint8Array(sig);
  const page = Deno.readTextFileSync("assets/pageops.html");
  Deno.writeFileSync(`${tss}-${a32h}.png`, fp_obj.im);
  Deno.writeTextFileSync(`${tss}-${a32h}.txt`, base64.bytesToBase64(u8sig));
  Deno.writeFileSync(`${backup}${tss}-${a32h}.png`, fp_obj.im);
  for await (const i of Deno.readDir("./")) {
    if (
      i.name != `${tss}-${a32h}.png` &&
      i.name.match(/^\d{8}T\d{9}Z\-\w{8}.png$/)
    ) {
      console.log(`removing ${i.name}`);
      Deno.remove(i.name);
    }
    if (
      i.name != `${tss}-${a32h}.txt` &&
      i.name.match(/^\d{8}T\d{9}Z\-\w{8}.txt$/)
    ) {
      console.log(`removing ${i.name}`);
      Deno.remove(i.name);
    }
  }
  Deno.writeTextFileSync(
    `${domain}.page.html`,
    page
      .replaceAll("thisisimage", `${tss}-${a32h}`)
      .replaceAll("thisisemoji", emoji)
      .replaceAll("thisislength",fp_obj.ln)
  );
}
function web_deal(req) {
  if (req.method == "GET") {
    const u = new URL(req.url);
    const page = u.pathname == "/"
      ? `triple.pub.page.html`
      : u.pathname.replace("/", "");
    let npg;
    let response;
    try {
      console.log(page);
      npg = Deno.readFileSync(page);
      const type = page.split(".").slice(-1);
      response = new Response(npg, {
        status: 200,
        headers: {
          "content-type": types[type],
        },
      });
    } catch {
      console.log("error 404");
      response = new Response(npg, {
        status: 404,
        headers: {
          "content-type": "text/plain;charset=utf-8",
        },
      });
    }
    return response;
  }
}
const types = {
  "js": "text/javascript;charset=utf-8",
  "css": "text/css",
  "svg": "image/svg+xml",
  "html": "text/html",
  "map": "application/json",
  "json": "application/json",
  "xz": "application/gzip",
  "png": "image/png",
  "zst": "application/zstd",
  "txt": "text/plain",
  "jpg": "image/jpg",
  "gif": "image/gif",
  "WebM": "video/webm",
  "mp4": "video/mp4",
  "mpg": "video/mp4",
  "webm": "video/webm",
  "ico": "image/x-icon",
};
Deno.serve({
  port: 3052,
  hostname: "0.0.0.0",
  handler: (req) => web_deal(req),
});
