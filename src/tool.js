import { create, web_deal } from "fpng-sign-serve";
const backup = Deno.env.get("CL_TRI_BACKUP");
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
create(site,backup)
Deno.serve({
  port: 3052,
  hostname: "0.0.0.0",
  handler: (req) => web_deal(req),
});